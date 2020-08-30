#' (Re)calculate rq for a given sample
#'
#' @param data A dataset output from pcr_tidy/pcr_control
#' @param relative_sample A sample to set others relative to (eg my_dmso_sample)
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
#' pcr_rq("U6D1")
#'
#' # Can also be run after using pcr_control:
#' pcr_tidy(dat_path) %>%
#' pcr_control("GAPDH") %>%
#' pcr_rq("U6D1")
pcr_rq <- function(data, relative_sample) {

        control_probe <- unique(data$control)

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
                              delta_ct_se  = .data$delta_ct_sd/sqrt(.data$rep),  # Rep might not be the correct metric here
                              df           = max(1, .data$rep + .data$rep[.data$target_name == control_probe] - 2),
                              t            = stats::qt(.05/2, .data$df, lower.tail = F)) %>%
                dplyr::group_by(.data$target_name) %>%
                dplyr::filter(!is.na(.data$sample_name)) %>%
                dplyr::mutate(delta_delta_ct = .data$delta_ct - .data$delta_ct[.data$sample_name == relative_sample],
                              rq     = 2^-(.data$delta_delta_ct),
                              rq_min = 2^-(.data$delta_delta_ct + .data$t*.data$delta_ct_se),
                              rq_max = 2^-(.data$delta_delta_ct - .data$t*.data$delta_ct_se)) %>%
                tidyr::unnest(.data$sample_nest)
}
