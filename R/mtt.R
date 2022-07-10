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
mtt <- function(x, ...) {
  UseMethod("mtt")
}

mtt.data.frame <- function(x, ...) {

}

mtt.gp <- function(x, ...) {

}

#' @export
#' @rdname mtt
mtt.spectramax <- function(spectramax, condition_names, drug_conc) {
  rs$data[[1]]$data |>
    gp::gp_sec("condition", nrow = 4, ncol = 6, labels = condition_names) |>
    gp::gp_sec("drug", nrow = 4, ncol = 1, labels = drug_conc, advance = FALSE) |>
    gp::gp_serve() |>
    dplyr::filter(!is.na(as.character(condition))) |>
    dplyr::group_by(.data$condition, .data$drug) |>
    dplyr::mutate(diff = .data$nm562 - .data$nm660,
                  mean = mean(.data$diff),
                  drug = as.numeric(levels(.data$drug)[.data$drug]),
                  drug = ifelse(.data$drug == 0, 1e-12, .data$drug)) |>
    dplyr::group_by(.data$condition) |>
    dplyr::mutate(div = .data$diff/.data$mean[which(.data$drug == min(.data$drug))])
}
