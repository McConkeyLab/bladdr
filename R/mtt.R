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
mtt_calc.data.frame <- function(x, ...) {
  x |>
    dplyr::filter(!is.na(as.character(condition))) |>
    dplyr::group_by(.data$condition, .data$drug) |>
    dplyr::mutate(diff = .data$nm562 - .data$nm660,
                  mean = mean(.data$diff),
                  drug = ifelse(.data$drug == 0, 1e-12, .data$drug)) |>
    dplyr::group_by(.data$condition) |>
    dplyr::mutate(div = .data$diff/.data$mean[which(.data$drug == min(.data$drug))])
}

mtt_calc.gp <- function(x, ...) {

}

#' @export
#' @rdname mtt_calc
mtt_calc.spectramax <- function(x, condition_names, drug_conc, ...) {
  x$data[[1]]$data |>
    gplate::gp_sec("condition", nrow = 4, ncol = 6, labels = condition_names) |>
    gplate::gp_sec("drug", nrow = 4, ncol = 1, labels = drug_conc, advance = FALSE) |>
    gplate::gp_serve() |>
    dplyr::filter(!is.na(as.character(condition))) |>
    dplyr::group_by(.data$condition, .data$drug) |>
    dplyr::mutate(diff = .data$nm562 - .data$nm660,
                  mean = mean(.data$diff),
                  drug = as.numeric(levels(.data$drug)[.data$drug]),
                  drug = ifelse(.data$drug == 0, 1e-12, .data$drug)) |>
    dplyr::group_by(.data$condition) |>
    dplyr::mutate(div = .data$diff/.data$mean[which(.data$drug == min(.data$drug))])
}



#' Plot MTT results
#'
#' @param mtt a `data.frame` output from  `mtt_calc()`
#'
#' @return a `ggplot`
#' @export
mtt_plot <- function(mtt) {
  mtt |>
    ggplot2::ggplot(ggplot2::aes(as.numeric(.data$drug),
                                 .data$div,
                                 color = .data$condition)) +
    ggplot2::geom_point() +
    ggplot2::scale_x_log10() +
    ggplot2::geom_smooth(method = drc::drm,
                         method.args = list(fct = drc::L.4(fixed = c(NA, NA, 1, NA))),
                         se = FALSE)
}
