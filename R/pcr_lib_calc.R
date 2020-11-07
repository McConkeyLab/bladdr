#' Calculate library PCR concentrations
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
#' pcr_lib_calc()

pcr_lib_calc <- function(tidy_pcr, dil_factor = 1000) {

        tidy_pcr %>%
                tidyr::nest(replicates = c(.data$well, .data$well_position, .data$ct, .data$quantity, .data$well_row, .data$well_col, dplyr::starts_with("prfdrop"), dplyr::starts_with("baxrox"))) %>%
                dplyr::group_by(.data$task) %>%
                dplyr::arrange(.data$ct_mean) %>%
                dplyr::mutate(standard_diff = .data$ct_mean - dplyr::lag(.data$ct_mean, default = .data$ct_mean[1]),
                              dil = 2^.data$standard_diff,
                              quant_actual = 6.8/cumprod(.data$dil),
                              dil = dplyr::if_else(.data$dil == 1, 0, .data$dil)) %>%
                tidyr::unnest(cols = .data$replicates) %>%
                dplyr::mutate(dil = dplyr::if_else(.data$task == "STANDARD", .data$dil, NA_real_),
                              standard_diff = dplyr::if_else(.data$task == "STANDARD", .data$standard_diff, NA_real_),
                              quant_actual = dplyr::if_else(.data$task == "STANDARD", .data$quant_actual, .data$quantity),
                              concentration = .data$quantity_mean * dil_factor)
}
