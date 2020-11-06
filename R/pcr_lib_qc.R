#' Generate visual library prep pcr quality control
#'
#' @param lib_calc_pcr an output from `pcr_lib_calc`
#'
#' @return a ggplot
#' @export
#'
#' @examples
#'
#' dat_path <- system.file("extdata", "untidy-standard-curve.xlsx", package = "bladdr")
#'
#' pcr_tidy(dat_path) %>%
#'         pcr_lib_calc() %>%
#'         pcr_lib_qc()
pcr_lib_qc <- function(lib_calc_pcr) {

        standards <- lib_calc_pcr %>%
                dplyr::filter(.data$task == "STANDARD")

        samples <- lib_calc_pcr %>%
                dplyr::filter(.data$task == "UNKNOWN")

        samples$quantity <- stats::predict(stats::lm(quantity~quant_actual, data = standards), samples)

        ggplot2::ggplot(standards, ggplot2::aes(x = .data$quant_actual, y = .data$quantity)) +
                ggplot2::geom_point(color = "#2297E6") +
                ggplot2::geom_smooth(method = "lm", color = "#2297E6", se = F) +
                ggplot2::geom_line(ggplot2::aes(x = .data$quantity, y = .data$quantity)) +
                ggplot2::geom_point(data = samples, ggplot2::aes(x = .data$quant_actual, y = .data$quantity), color = "#61D04F", alpha = 0.7) +
                ggplot2::geom_point(ggplot2::aes(x = .data$quantity, y = .data$quantity)) +
                ggplot2::scale_x_log10() +
                ggplot2::scale_y_log10()
}
