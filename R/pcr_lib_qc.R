#' Title
#'
#' @param libcalc_pcr
#'
#' @return
#' @export
#'
#' @examples
#'
dat_path <- system.file("extdata", "untidy-standard-curve.xlsx", package = "bladdr")

pcr_tidy(dat_path) %>%
        pcr_lib_calc() %>%
        pcr_lib_qc()
pcr_lib_qc <- function(lib_calc_pcr) {

        standards <- lib_calc_pcr$standards %>%
                dplyr::mutate(
                        dil_factor = dplyr::case_when(
                                is.na(.data$dil) ~ 1,
                                TRUE ~ .data$dil),
                        dil_running = 6.8/cumprod(.data$dil_factor),
                        dil_theoretical = 6.8 * 10^-(0:4)
                )
        print(standards)
        ggplot2::ggplot(standards, ggplot2::aes(x = dil_running, y = dil_theoretical)) +
                ggplot2::geom_point(data = lib_calc_pcr$samples, aes = c(x = quantity_mean, y = quantity_mean)) +
                ggplot2::geom_point() +
                ggplot2::geom_smooth(method = "lm", se = F) +
                ggplot2::geom_line(aes(x = dil_theoretical, y = dil_theoretical)) +
                ggplot2::scale_x_log10() +
                ggplot2::scale_y_log10()
}
