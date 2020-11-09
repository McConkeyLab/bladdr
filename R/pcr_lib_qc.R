#' Generate visual library prep pcr quality control
#'
#' @param lib_calc_pcr an output from `pcr_lib_calc`
#'
#' @return a ggplot
#' @export
#'
#' @importFrom ggplot2 aes element_blank
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

        sample_summary <- samples %>%
                dplyr::group_by(.data$sample_name) %>%
                dplyr::summarize(quantity_mean = mean(.data$quantity_mean))

        standard_summary <- standards %>%
                dplyr::group_by(.data$quantity) %>%
                dplyr::summarize(quantity_mean = mean(.data$quantity),
                                 quant_actual = mean(.data$quant_actual),
                                 dil = mean(.data$dil)) %>%
                tidyr::pivot_longer(cols = c(.data$quantity_mean, .data$quant_actual))

        dilution_lines <- standard_summary %>%
                dplyr::filter(.data$name == "quant_actual") %>%
                dplyr::mutate(line_start = 1/.data$value,
                              line_end = dplyr::lag(.data$line_start),
                              dil = dplyr::lag(.data$dil),
                              y = rep_len(c(1.1, 0.9), 5),
                              y_text = rep_len(c(1.15, 0.85), 5)) %>%
                dplyr::filter(!is.na(.data$line_end)) %>%
                dplyr::rowwise() %>%
                dplyr::mutate(mid = sqrt(.data$line_start*.data$line_end))

        vert_lines <- tibble::tibble(x = c(dilution_lines$line_start,
                                           dilution_lines$line_end)) %>%
                dplyr::arrange(.data$x) %>%
                dplyr::mutate(y = rep(c(1.1, 0.9, 1.1, 0.9), each = 2),
                              yend = rep(c(1.05, 0.95, 1.05, 0.95), each = 2))

        ggplot2::ggplot(standard_summary, aes(x = 1/.data$value, y = 1, color = .data$name)) +
                ggplot2::geom_point(size = 10, alpha = 0.7) +
                ggplot2::scale_color_manual(values = c("#00AAAA", "#222222")) +
                ggplot2::geom_segment(data = dilution_lines, aes(x = .data$line_start, xend = .data$line_end, y = .data$y, yend = .data$y), size = 1) +
                ggplot2::geom_segment(data = vert_lines, aes(x = .data$x, xend = .data$x, y = .data$y, yend = .data$yend), color = "#00AAAA", size = 1) +
                ggplot2::geom_text(data = dilution_lines, aes(label = round(.data$dil, 1), y = .data$y_text, x = .data$mid), size = 8) +
                ggplot2::scale_x_log10() +
                ggplot2::theme_minimal() +
                ggplot2::theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "none") +
                ggplot2::coord_cartesian(ylim = c(0.5, 1.5)) +
                ggplot2::geom_point(data = sample_summary, aes(x = 1/.data$quantity_mean, y = 1), color = "red")
}
