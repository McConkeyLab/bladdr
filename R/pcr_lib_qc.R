#' Generate visual library prep pcr quality control report
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
#' lib_calc_pcr <- pcr_tidy(dat_path, pad_zero = TRUE) %>%
#'         pcr_lib_calc() %>%
#'         pcr_lib_qc()
#'
pcr_lib_qc <- function(lib_calc_pcr) {
        # Does mean(quantity) == quantity_mean?
        dat <- lib_calc_pcr %>%
                dplyr::select(.data$task, .data$sample_name, .data$quantity_mean, .data$concentration, .data$quantity, .data$quant_actual, .data$dil, .data$slope, .data$efficiency, .data$r2, .data$ct)

        standards <- dat %>%
                dplyr::filter(.data$task == "STANDARD")

        samples <- dat %>%
                dplyr::filter(.data$task == "UNKNOWN")

        sample_summary <- samples %>%
                dplyr::group_by(.data$sample_name) %>%
                dplyr::summarize(quantity_mean = mean(.data$quantity_mean),
                                 concentration_mean = mean(.data$concentration))

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

        standards_plot <-
                ggplot2::ggplot(standard_summary, aes(x = 1/.data$value, y = 1, color = .data$name)) +
                ggplot2::geom_point(size = 10, alpha = 0.7) +
                ggplot2::scale_color_manual(values = c("#00AAAA", "#222222")) +
                ggplot2::geom_segment(data = dilution_lines, aes(x = .data$line_start, xend = .data$line_end, y = .data$y, yend = .data$y), size = 1) +
                ggplot2::geom_segment(data = vert_lines, aes(x = .data$x, xend = .data$x, y = .data$y, yend = .data$yend), color = "#00AAAA", size = 1) +
                ggplot2::geom_text(data = dilution_lines, aes(label = round(.data$dil, 1), y = .data$y_text, x = .data$mid), size = 8) +
                ggplot2::scale_x_log10() +
                ggplot2::theme_minimal() +
                ggplot2::theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "none") +
                ggplot2::coord_cartesian(ylim = c(0.8, 1.2)) +
                ggplot2::geom_point(data = sample_summary, aes(x = 1/.data$quantity_mean, y = 1), color = "red")

        slope_text <-
                dplyr::tibble(x = -2.5,
                              y = 20,
                              label = paste(paste("Slope:", round(standards$slope[1], 2)),
                                            paste("R\u00B2:", round(standards$r2[1], 2)),
                                            paste0("Efficiency: ", round(standards$efficiency[1], 1), "%"),
                                            sep = "\n"))

        slope_plot <-
                ggplot2::ggplot(standards, aes(x = log10(.data$quantity), y = .data$ct)) +
                ggplot2::geom_point(color = "#00AAAA", size = 5) +
                ggplot2::geom_smooth(method = "lm", se = F, color = "#666666") +
                ggplot2::geom_text(data = slope_text, aes(x = .data$x,  y = .data$y, label = .data$label), size = 8) +
                ggplot2::xlab("Log10(Quantity)") +
                ggplot2::ylab("CT")

        outliers_plot_dat <- dat %>%
                tidyr::nest(reps = c(.data$ct, .data$quantity, .data$quant_actual, .data$concentration)) %>%
                dplyr::mutate(mean_wo_outlier = purrr::map(.data$reps, find_mean)) %>%
                tidyr::unnest_wider(.data$mean_wo_outlier) %>%
                tidyr::unnest(cols = c(.data$no_po_mean, .data$no_po_sd, .data$reps)) %>%
                dplyr::mutate(keep = dplyr::case_when(.data$no_po_mean - (3*.data$no_po_sd) < .data$ct & .data$no_po_mean + (3*.data$no_po_sd) > .data$ct ~ T,
                                                      is.na(.data$no_po_sd)~T,
                                                      T~NA),
                              sample_name = dplyr::case_when(is.na(.data$sample_name) & .data$quantity > 6.0000 ~ "1",
                                                             is.na(.data$sample_name) & .data$quantity > 0.6000 ~ "1:10",
                                                             is.na(.data$sample_name) & .data$quantity > 0.0600 ~ "1:100",
                                                             is.na(.data$sample_name) & .data$quantity > 0.0060 ~ "1:1000",
                                                             is.na(.data$sample_name) & .data$quantity > 0.0006 ~ "1:10000",
                                                             T ~   .data$sample_name)) %>%
                dplyr::group_by(.data$sample_name) %>%
                dplyr::mutate(adj_mean = mean(.data$keep*.data$ct, na.rm = T),
                              adj_sd   = stats::sd(.data$keep*.data$ct, na.rm = T),
                              keep_logi = dplyr::if_else(!is.na(.data$keep), T, F),
                              z = (.data$ct-.data$adj_mean)/.data$adj_sd,
                              overflow = abs(.data$z) > 10,
                              z_plot = dplyr::case_when(.data$z > 10 ~ 10,
                                                        .data$z < -10 ~ -10,
                                                        T ~ .data$z),
                              label = dplyr::case_when(.data$z > 10 ~ paste(as.character(round(.data$z, 0)),">>>"),
                                                       .data$z < -10 ~ paste("<<<", as.character(round(.data$z, 0))),
                                                       T ~ NA_character_)) %>%
                dplyr::filter(!is.na(.data$sample_name))

        outliers_plot <-
                ggplot2::ggplot(outliers_plot_dat, aes(x = .data$z_plot, y = .data$sample_name, color = .data$keep_logi, shape = .data$overflow)) +
                ggplot2::scale_color_manual(values = c("#666666", "#00AAAA")) +
                ggplot2::scale_shape_manual(values = c(16, NA)) +
                ggplot2::geom_text(aes(label = .data$label), size = 5) +
                ggplot2::geom_rect(xmin = -3, xmax = 3, ymin = 0, ymax = nrow(outliers_plot_dat), fill = "#AACCCC",  alpha = 0.1,  color = NA) +
                ggplot2::geom_vline(xintercept = 0, color = "#555555") +
                ggplot2::geom_point(size = 5) +
                ggplot2::xlab("Z-score") +
                ggplot2::coord_cartesian(xlim = c(-11.5, 11.5), ylim = c(0.5, length(unique(outliers_plot_dat$sample_name)) + 0.5), expand = F) +
                ggplot2::theme(panel.grid.major.x = element_blank(), axis.title.y = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none")


        table_samples <- sample_summary %>%
                dplyr::select(-.data$quantity_mean) %>%
                dplyr::rename("Sample Name" = .data$sample_name,
                              "Concentration" = .data$concentration_mean)

        conc_plot <-
                ggplot2::ggplot(table_samples, aes(x = .data$`Sample Name`, y = .data$`Concentration`)) +
                ggplot2::geom_point(color = "#00AAAA", size = 5) +
                ggplot2::theme(axis.title = element_blank(), legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))

        report <- system.file("lib-qc.Rmd", package = "bladdr")
        rmarkdown::render(report, params = list(stan = standards_plot, out = outliers_plot, samp = table_samples, conc = conc_plot, slope = slope_plot))
}


find_mean <- function(df){
        if(nrow(df) >= 3 & !all(is.na(df$ct))) {
                hc <- df$ct %>%
                        stats::dist() %>%
                        stats::hclust()
                possible_outlier <- hc$merge[nrow(df) - 1, 1] %>%
                        abs()
                no_po <- df$ct[-possible_outlier]
                no_po_mean <- mean(no_po)
                no_po_sd <- stats::sd(no_po)
                list(no_po_mean = no_po_mean, no_po_sd = no_po_sd)
        } else {
                no_po_mean <- mean(df$ct)
                no_po_sd <- stats::sd(df$ct, na.rm = T)
                list(no_po_mean = no_po_mean, no_po_sd = no_po_sd)
        }
}
