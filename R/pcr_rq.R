#' (Re)calculate rq for a given sample
#'
#' @param data A dataset output from pcr_tidy/pcr_control
#' @param relative_sample A sample to set others relative to (eg my_dmso_sample)
#'
#' @return A tibble
#' @export
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
                dplyr::distinct(`Sample Name`, `Target Name`, .keep_all = T) %>%
                dplyr::filter(!is.na(`Sample Name`)) %>%
                dplyr::group_by(`Target Name`) %>%
                dplyr::mutate(d_ct_mean = as.numeric(`Delta Ct Mean`)) %>%
                dplyr::mutate(ddct = d_ct_mean - d_ct_mean[`Sample Name` == relative_sample],
                              rerq = 1/2^ddct)
}
