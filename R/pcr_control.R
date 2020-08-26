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
                dplyr::distinct(.data$`Sample Name`, .data$`Target Name`, .keep_all = T) %>%
                dplyr::filter(!is.na(.data$`Sample Name`)) %>%
                dplyr::group_by(.data$`Sample Name`) %>%
                dplyr::mutate(`Delta Ct Mean` = .data$`Ct Mean` - .data$`Ct Mean`[.data$`Target Name` == control_probe])
}
