#' (Re)calculate Delta Ct Mean
#'
#' @param data A dataset output from pcr_tidy
#' @param control_probe A probe to be used as an endogenous control (eg GAPDH)
#'
#' @return A tibble
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr")
#'
#' pcr_tidy(dat_path) %>%
#' pcr_control("GAPDH")
pcr_control <- function(data, control_probe) {
        data %>%
                dplyr::group_by(.data$target_name, .data$sample_name) %>%
                dplyr::mutate(ct_mean = mean(.data$ct),
                              ct_sd   = stats::sd(.data$ct),
                              rep     = dplyr::n()) %>%
                dplyr::ungroup() %>%
                tidyr::nest(sample_nest = c(.data$ct, .data$well, .data$well_row, .data$well_col, .data$baseline_start, .data$baseline_end)) %>%
                dplyr::group_by(.data$sample_name) %>%
                dplyr::mutate(delta_ct     = .data$ct_mean - .data$ct_mean[.data$target_name == control_probe],
                              delta_ct_sd  = sqrt(.data$ct_sd^2 + .data$ct_sd[.data$target_name == control_probe]^2),
                              delta_ct_se    = .data$delta_ct_sd/sqrt(.data$rep),  # Rep might not be the correct metric here
                              df           = max(1, .data$rep + .data$rep[.data$target_name == control_probe] - 2),
                              t            = stats::qt(.05/2, .data$df, lower.tail = F)) %>%
                tidyr::unnest(.data$sample_nest)
}
