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
#' a <- pcr_tidy(dat_path) %>%
#' pcr_lib_calc()
#'
pcr_lib_calc <- function(tidy_pcr, dil_factor = 1000) {

        standards <- tidy_pcr %>%
                dplyr::filter(.data$task == "STANDARD") %>%
                tidyr::nest(replicates = c(well, well_position, ct, well_row, well_col)) %>%
                dplyr::group_by(task) %>%
                dplyr::mutate(standard_diff = ct_mean - dplyr::lag(ct_mean, default = ct_mean[1]))
                # dplyr::group_by(.data$quantity) %>%
                # dplyr::arrange(.data$ct_mean) %>%
                # dplyr::group_by(.data$quantity, .data$ct_mean) %>%
                # dplyr::mutate(standard_diff = .data$ct_mean - dplyr::lag(.data$ct_mean))
                # dplyr::mutate(dil = 2^.data$standard_diff)

        # samples <- tidy_pcr %>%
        #         dplyr::filter(.data$task == "UNKNOWN") %>%
        #         dplyr::group_by(.data$sample_name) %>%
        #         dplyr::summarize(ct_mean = mean(.data$ct),
        #                          quantity_mean = mean(.data$quantity),
        #                          concentration = .data$quantity_mean * dil_factor)
        #
        # output <- list(input = tidy_pcr, standards = standards, samples = samples)
}
