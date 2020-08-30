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
                dplyr::distinct(.data$sample_name, .data$target_name, .keep_all = T) %>%
                dplyr::filter(!is.na(.data$sample_name)) %>%
                dplyr::group_by(.data$sample_name) %>%
                dplyr::mutate(delta_ct_mean = .data$ct_mean - .data$ct_mean[.data$target_name == control_probe])
}
