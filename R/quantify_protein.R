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
               remove_empty = TRUE,
               remove_outliers = c("all", "samples", "standards", "none")) {

  replicate_orientation <- rlang::arg_match(replicate_orientation)
  remove_outliers <- rlang::arg_match(remove_outliers)

  abs <- qp_read(x)
  abs_tidy <- qp_tidy(abs, replicate_orientation)
  mean_abs <- qp_calc_abs_mean(abs_tidy, remove_outliers)
  fit <- qp_fit(mean_abs)
  conc <- qp_calc_conc(mean_abs, fit)

  if (remove_empty) {
    conc <- dplyr::filter(
      conc, .data$.pred_conc > 0 | .data$sample_type == "standard"
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

qp_calc_abs_mean <- function(x, remove_outliers) {
  standards <- x |>
    dplyr::filter(.data$sample_type == "standard") |>
    calc_mean(remove_outliers %in% c("all", "standards"))
  unknowns <- x |>
    dplyr::filter(.data$sample_type == "unknown") |>
    calc_mean(remove_outliers %in% c("all", "samples"))
  rbind(standards, unknowns)
}

calc_mean <- function(df, remove_outliers) {
  df <- dplyr::group_by(df, .data$sample_type, .data$index)
  if (remove_outliers) {
    df <- df |>
      dplyr::mutate(
        is_outlier = mark_outlier(.data$value),
        mean = mean(.data$value[!.data$is_outlier], na.rm = TRUE)
      )
  } else {
    df <- df |>
      dplyr::mutate(
        is_outlier = NA,
        mean = mean(.data$value, na.rm = TRUE)
      )
  }
  #FIXME log calculation probably does not belong here
  df <- df |>
    dplyr::mutate(log_abs = log2(.data$value)) |>
    dplyr::ungroup()
}

mark_suspect <- function(nums) {
  # Marking a suspect with 2 or fewer samples doesn't make sense
  na_index <- which(!is.na(nums))
  no_na <- na.omit(nums)
  if (length(no_na) <= 2) return(rep(FALSE, length(nums)))
  hc <- stats::hclust(stats::dist(no_na))
  no_na_index <- abs(hc$merge[length(no_na) -1, 1])
  suspect_index <- na_index[no_na_index]
  out <- rep(FALSE, length(nums))
  out[suspect_index] <- TRUE
  out
}

mark_outlier <- function(nums) {
  marked <- mark_suspect(nums)
  if (!any(marked)) return(marked)
  no_suspect <- nums[!marked]
  suspect <- nums[marked]
  mean_no_suspect <- mean(no_suspect, na.rm = TRUE)
  sd_no_suspect <- sd(no_suspect, na.rm = TRUE)
  suspect_is_outlier <- abs(suspect - mean_no_suspect) > (3 * sd_no_suspect)
  if (suspect_is_outlier) {
    return(marked)
  } else {
    return(rep(FALSE, length(nums)))
  }
}

# Fit conc ~ abs using standards absorbances -----------------------------------
qp_fit <- function(x) {
  standards <- x |>
    dplyr::filter(
      .data$sample_type == "standard",
      f_or_na(.data$is_outlier)
    )
  fit_data <- dplyr::mutate(standards, log_conc = log2(.data$conc + .5))
  stats::lm(log_conc ~ log_abs, data = fit_data)
}

# Predict concentrations from standards fit ------------------------------------
qp_calc_conc <- function(x, fit) {
  with_predictions <- dplyr::bind_cols(x, .pred = stats::predict(fit, x))
  with_predictions |>
    dplyr::mutate(.pred_conc = (2^.data$.pred) - 0.5) |>
    dplyr::group_by(.data$sample_type, .data$index) |>
    dplyr::mutate(
      .pred_conc_mean = mean(.data$.pred_conc[which(f_or_na(.data$is_outlier))])
    ) |>
    dplyr::ungroup()
}

f_or_na <- function(x) {
  !x | is.na(x)
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

  plot_data <- x$qp |>
    dplyr::mutate(
      outlier = !f_or_na(.data$is_outlier),
      y = ifelse(
        .data$sample_type == "standard",
        log2(.data$conc + 0.5),
        log2(.data$.pred_conc_mean + 0.5)
      )
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$log_abs, y = .data$y,
      color = .data$sample_type,
      shape = .data$outlier
    )) +
    ggplot2::scale_color_viridis_d(
               option = "viridis", end = 0.8, direction = -1
             ) +
    ggplot2::geom_abline(intercept = x$fit$coefficients[1],
                slope = x$fit$coefficients[2],
                size = 2,
                alpha = 0.2) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::scale_shape_manual(values = c(16, 4)) +
    ggplot2::labs(x = "Log2(Absorbance)", y = "Log2(Concentration + 0.5)")
}
