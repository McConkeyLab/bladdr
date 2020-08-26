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
        data %>%
                dplyr::distinct(.data$`Sample Name`, .data$`Target Name`, .keep_all = T) %>%
                dplyr::filter(!is.na(.data$`Sample Name`)) %>%
                dplyr::group_by(.data$`Target Name`) %>%
                dplyr::mutate(`Delta Delta Ct` = .data$`Delta Ct Mean` - .data$`Delta Ct Mean`[.data$`Sample Name` == relative_sample],
                              RQ = 1/2^.data$`Delta Delta Ct`)
}
