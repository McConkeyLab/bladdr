#' (Re)calculate Delta Ct Mean
#'
#' @param data A dataset output from pcr_tidy
#' @param control_probe A probe to be used as an endogenous control (eg GAPDH)
#'
#' @return A tibble
#' @export
#'
#' @examples
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr")
#'
#' pcr_tidy(dat_path) %>%
#' pcr_control("GAPDH")
pcr_control <- function(data, control_probe) {
        data %>%
                dplyr::distinct(`Sample Name`, `Target Name`, .keep_all = T) %>%
                dplyr::filter(!is.na(`Sample Name`)) %>%
                dplyr::group_by(`Sample Name`) %>%
                dplyr::mutate(`Delta Ct Mean` = `Ct Mean` - `Ct Mean`[`Target Name` == control_probe])
}
