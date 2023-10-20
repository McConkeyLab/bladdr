#' Quantify protein concentration from a MicroBCA assay
#'
#' @param x A `spectramax`, `gp`, or `data.frame` object, or path to SPECTRAmax
#' .xls(x)/.txt file.
#' @param replicate_orientation Either 'h' or 'v' - see Details.
#' @param sample_names Character vector of sample names.
#' @param remove_empty Should wells that have less absorbance than the lowest
#' standard be dropped?
#'
#' @details The standards must be in ascending concentration starting in the
#'   upper left corner. Whether this is from from left to right or top to bottom
#'   can be specified in 'replicate orientation'. Note that 'replicate
#'   orientation' specified the direction that REPLICATES lie, NOT the direction
#'   the samples flow (which will be opposite).
#'
#' @return a `tibble`
#' @export
#' @importFrom rlang .data
qp <- function(x,
               replicate_orientation = c("h", "v"),
               sample_names = NULL,
               remove_empty = TRUE) {

  replicate_orientation <- rlang::arg_match(replicate_orientation)

  abs <- qp_read(x)

  # Derived constants
  max_samples <- gplate::wells(abs) %/% 3
  n_standards <- 7
  max_unknowns <- max_samples - n_standards

  mean_abs <- abs |>
    qp_tidy(replicate_orientation, max_unknowns) |>
    qp_calc_abs_mean()

  fit <- mean_abs |>
    dplyr::filter(.data$sample_type == "standard") |>
    dplyr::filter(.data$keep) |>
    qp_fit()

  conc <- qp_calc_conc(mean_abs, fit)

  if (remove_empty) {
    conc <-
      dplyr::filter(
               conc, .data$pred_conc > 0 | .data$sample_type == "standard"
             )
  }

  if (!is.null(sample_names)) {
    # Will return "NA" instead of erroring if sample names < # samples
    length(sample_names) <- max(conc$index, na.rm = TRUE)
  } else {
    sample_names <- as.character(1:max(conc$index, na.rm = TRUE))
  }

  qp <- conc |>
    dplyr::mutate(sample_name = ifelse(
                    .data$sample_type == "unknown",
                    sample_names[.data$index],
                    paste("Standard", .data$index))
                  )

  list(fit = fit, qp = qp, gp = abs)
}

# Tidy qp gp according to replicate orientation --------------------------------
qp_tidy <- function(x, replicate_orientation) {
  max_samples <- gplate::wells(x) %/% 3
  n_standards <- 7
  max_unknowns <- max_samples - n_standards

  if (replicate_orientation == "v") {
    nrow <- nrow2 <- 3
    ncol <- c(7, max_unknowns)
    flow <- "row"
    ncol2 <- 1
  } else {
    nrow <- c(7, max_unknowns)
    ncol <- ncol2 <- 3
    flow <- "col"
    nrow2 <- 1
  }

  x |>
    gplate::gp_sec(name = "sample_type", nrow, ncol, wrap = TRUE, flow = flow,
               labels = c("standard", "unknown"), break_sections = FALSE) |>
    gplate::gp_sec(name = "index", nrow2, ncol2, break_sections = FALSE) |>
    gplate::gp_serve() |>
    dplyr::mutate(
             index = as.numeric(.data$index),
             conc = ifelse(.data$index > 1, 2^(.data$index - 5), 0),
             conc = ifelse(
               .data$sample_type == "standard",
               .data$conc,
               NA_real_
             )
           )
}


# Calculate outlier-free absorbance means --------------------------------------

qp_calc_abs_mean <- function(x) {

  x |>
    dplyr::group_by(.data$sample_type, .data$index, .data$conc) |>
    tidyr::nest() |>
    dplyr::mutate(mean_no_outlier = purrr::map(.data$data, find_mean)) |>
    dplyr::select(-.data$data) |>
    tidyr::unnest(.data$mean_no_outlier) |>
    dplyr::group_by(.data$sample_type, .data$index) |>
    dplyr::mutate(
             no_out_mean = mean(
               .data$value[!.data$is_suspect],
               na.rm = TRUE
             ),
             no_out_sd = stats::sd(
                                  .data$value[!.data$is_suspect],
                                  na.rm = TRUE
                                ),
             ##FIXME This SHOULD work if I just ! the output of the first arg
             keep = ifelse(
               abs(.data$value - .data$no_out_mean) > 3 * .data$no_out_sd,
               FALSE,
               TRUE
             ),
                  log_abs = log2(.data$value)) |>
    dplyr::ungroup()
}

find_mean <- function(df) {
  sample_hclust <- df$value |>
    stats::dist() |>
    stats::hclust()
  suspect_index <- sample_hclust$merge[nrow(df) - 1, 1] |>
    abs()
  df$is_suspect <- FALSE
  df$is_suspect[suspect_index] <- TRUE
  tibble::tibble(value = df$value, is_suspect = df$is_suspect)
}

# Fit conc ~ abs using standards absorbances -----------------------------------
qp_fit <- function(standards) {
  fit_data <- standards |>
    dplyr::mutate(log_conc = log2(.data$conc + .5))

  stats::lm(log_conc ~ log_abs, data = fit_data)
}

# Predict concentrations from standards fit ------------------------------------
qp_calc_conc <- function(x, fit) {
  x |>
    dplyr::bind_cols(.pred = stats::predict(fit, x)) |>
    tidyr::unnest(dplyr::everything()) |>
    dplyr::group_by(.data$sample_type, .data$index) |>
    dplyr::mutate(true_mean = ifelse(
                    .data$sample_type == "standard",
                    log2(.data$conc + 0.5),
                    mean(.data$.pred[.data$keep])
                  ),
                  pred_conc = (2^.data$true_mean) - .5) |>
    dplyr::ungroup()
}

# Calculate Dilutions ----------------------------------------------------------
#' Calculate dilutions for an analyzed `qp` `list`
#'
#' @param x The output of `qp()`
#' @param target_conc Target concentration in (mg/mL) protein
#' @param target_vol Target volume in uL
#'
#' @return A list, where the `qp` item has volumes of lysate and volumes of H2O
#' to add.
#' @export
qp_calc_dil <- function(x, target_conc, target_vol) {
  x$qp <- x$qp |>
    dplyr::rowwise() |>
    dplyr::mutate(
             .temp = list(
               dilute(.data$pred_conc, target_conc, target_vol, quiet = TRUE)
             )
           ) |>
    tidyr::unnest_wider(.data$.temp)
  x
}

# Visualization ----------------------------------------------------------------
#' View the absorbances of an analyzed `qp` as they were on the plate
#'
#' @param x The output of `qp()`
#'
#' @return a `ggplot`
#' @export
make_qp_plate_view <- function(x) {
  x$gp |>
    gplate::gp_plot(.data$value) +
    ggplot2::geom_point(ggplot2::aes(color = .data$value), size = 20) +
    ggplot2::geom_text(
               ggplot2::aes(label = round(.data$value, 2)),
               color = "black"
             ) +
    ggplot2::scale_color_gradient(low = "darkseagreen1", high = "mediumpurple3")
}

#' View an absorbance/concentration plot
#'
#' @param x The output of `qp()`
#'
#' @return a `ggplot`
#' @export
make_qp_standard_plot <- function(x) {
  ggplot2::ggplot(x$qp, ggplot2::aes(x = .data$log_abs,
                   y = .data$true_mean,
                   color = .data$sample_type,
                   shape = .data$keep)) +
    ggplot2::scale_color_viridis_d(
               option = "viridis", end = 0.8, direction = -1
             ) +
    ggplot2::geom_abline(intercept = x$fit$coefficients[1],
                slope = x$fit$coefficients[2],
                size = 2,
                alpha = 0.2) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::scale_shape_manual(values = c(4, 16)) +
    ggplot2::labs(x = "Log2(Absorbance)", y = "Log2(Concentration + 0.5)")
}
