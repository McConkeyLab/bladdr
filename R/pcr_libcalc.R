#' Calculate library pcr concentrations
#'
#' @param tidy_pcr a dataset run that has been run through `pcr_tidy()`
#' @param dil_factor the factor to which the libraries were diluted for pcr
#'
#' @return a list, containing the original dataframe, data on the standards, and data on the samples.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#'
#' dat_path <- system.file("extdata", "untidy-standard-curve.xlsx", package = "bladdr")
#'
#' pcr_tidy(dat_path) %>%
#' pcr_libcalc()
#'
pcr_libcalc <- function(tidy_pcr, dil_factor = 1000) {

        standards <- tidy_pcr %>%
                dplyr::filter(.data$task == "STANDARD") %>%
                dplyr::group_by(.data$quantity) %>%
                dplyr::summarize(ct_mean = mean(.data$ct),
                                 slope = unique(.data$slope),
                                 efficiency = unique(.data$efficiency)) %>%
                dplyr::arrange(.data$ct_mean) %>%
                dplyr::mutate(standard_diff = .data$ct_mean - dplyr::lag(.data$ct_mean)) %>%
                dplyr::mutate(dil = 2^.data$standard_diff)

        samples <- tidy_pcr %>%
                dplyr::filter(.data$task == "UNKNOWN") %>%
                dplyr::group_by(.data$sample_name) %>%
                dplyr::summarize(ct_mean = mean(.data$ct),
                                 quantity_mean = mean(.data$quantity),
                                 concentration = .data$quantity_mean * dil_factor)

        output <- list(tidy_pcr, standards, samples)
}
