#' Make a plot from MTT data
#'
#' @param x A `data.frame`, `gp`, or `spectramax` object. See Details.
#' @param ... Arguments passed to their repsective methods
#'
#' @details If supplied with a `data.frame`, `mtt` will expect...
#'
#' Using a `gp` object is a good idea if you have a 'non-standard' plate layout
#' (standard being each quarter of the plate is a condition)
#'
#' @return a `ggplot`
#' @export
mtt_calc <- function(x, ...) {
  UseMethod("mtt_calc")
}

# Assume there are already condition and drug columns
#' @export
#' @rdname mtt_calc
mtt_calc.data.frame <- function(x, drug_conc, ic_pct = 50, ...) {
  x |>
    dplyr::filter(!is.na(as.character(.data$condition))) |>
    dplyr::group_by(.data$condition, .data$drug) |>
    dplyr::mutate(
      diff = .data$nm562 - .data$nm660,
      mean = mean(.data$diff),
      drug = ifelse(.data$drug == 0, 1e-12, .data$drug)
    ) |>
    dplyr::group_by(.data$condition) |>
    dplyr::mutate(div = .data$diff / .data$mean[which(.data$drug == min(.data$drug))]) |>
    dplyr::group_by(.data$condition) |>
    tidyr::nest() |>
    dplyr::mutate(
      fit = purrr::map(.data$data, mtt_model),
      curve = purrr::map(.data$fit, mtt_make_curve, drug_conc),
      ic = purrr::map(.data$fit, mtt_get_ic, ic_pct = ic_pct)
    ) |>
    tidyr::unnest(cols = c(.data$ic, .data$data)) |>
    dplyr::ungroup()
}

mtt_calc.gp <- function(x, ...) {

}

#' @export
#' @rdname mtt_calc
mtt_calc.spectramax <- function(x, condition_names, drug_conc, ic_pct = 50, ...) {
  x$data[[1]]$data |>
    gplate::gp_sec("condition", nrow = 4, ncol = 6, labels = condition_names) |>
    gplate::gp_sec("drug", nrow = 4, ncol = 1, labels = drug_conc, advance = FALSE) |>
    gplate::gp_serve() |>
    dplyr::filter(!is.na(as.character(.data$condition))) |>
    dplyr::group_by(.data$condition, .data$drug) |>
    dplyr::mutate(
      diff = .data$nm562 - .data$nm660,
      mean = mean(.data$diff),
      drug = as.numeric(levels(.data$drug)[.data$drug]),
      drug = ifelse(.data$drug == 0, 0, .data$drug)
    ) |>
    dplyr::group_by(.data$condition) |>
    dplyr::mutate(div = .data$diff / .data$mean[which(.data$drug == min(.data$drug))]) |>
    dplyr::group_by(.data$condition) |>
    tidyr::nest() |>
    dplyr::mutate(
      fit = purrr::map(.data$data, mtt_model),
      curve = purrr::map(.data$fit, mtt_make_curve, drug_conc),
      ic = purrr::map(.data$fit, mtt_get_ic, ic_pct = ic_pct)
    ) |>
    tidyr::unnest(cols = c(.data$ic, .data$data)) |>
    dplyr::ungroup()
}

mtt_model <- function(data) {
  drc::drm(
    div ~ drug,
    data = data, fct = drc::LL.4(fixed = c(NA, NA, NA, NA)),
    lowerl = c(-Inf, 0, -Inf, -Inf)
  )
}

mtt_make_curve <- function(fit, drug_conc, length_out = 1000) {
  min <- ifelse(min(drug_conc) == 0, 1e-12, min(drug_conc))
  x <- exp(seq(log(min), log(max(drug_conc)), length.out = length_out))
  curve <- fit$curve[[1]](x)
  data.frame(x = x, y = curve)
}

mtt_get_ic <- function(fit, ic_pct) {
  drc::ED(fit, respLev = ic_pct, display = FALSE) |> #
    dplyr::as_tibble() |>
    dplyr::rename(
      ic_value = .data$Estimate,
      ic_std_err = .data$`Std. Error`
    ) |>
    dplyr::mutate(ic_pct = ic_pct)
}

#' Plot MTT results
#'
#' @param mtt a `data.frame` output from  `mtt_calc()`
#' @param plot_ics logical. Should the calculated inhibitor concentrations be
#'   plotted?
#'
#' @return a `ggplot`
#' @export
mtt_plot <- function(mtt, plot_ics = FALSE) {
  fit_curve <- mtt |>
    dplyr::select("condition", "curve") |>
    tidyr::unnest(.data$curve)

  plot <- mtt |>
    ggplot2::ggplot(
      ggplot2::aes(
        as.numeric(.data$drug),
        .data$div,
        color = .data$condition
      )
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_x_log10() +
    ggplot2::geom_line(data = fit_curve, ggplot2::aes(.data$x, .data$y))

  if (plot_ics) {
    ic_annot <- mtt |>
      dplyr::group_by(.data$condition) |>
      dplyr::summarize(
        ic_value = unique(.data$ic_value),
        ic_std_err = unique(.data$ic_std_err),
        ic_pct = unique(.data$ic_pct),
        ic_fit = unique(.data$fit)
      ) |>
      dplyr::mutate(
        label = paste0("IC", ic_pct, ": ", round(ic_value, 2))
      ) |>
      dplyr::mutate(
        ic_y_fn = purrr::map(.data$ic_fit, \(x) x[[3]][[1]])
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(ic_y = .data$ic_y_fn(.data$ic_value)[, 1]) |>
      dplyr::ungroup()
    plot <- plot +
      ggrepel::geom_label_repel(
        data = ic_annot,
        ggplot2::aes(x = .data$ic_value, y = .data$ic_y, label = .data$label),
        nudge_x = -5
      )
  }
  plot
}
