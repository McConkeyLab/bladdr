#' Plot qPCR results
#'
#' @param tidy_pcr an output from the `pcr_tidy` function, or some derivative thereof
#'
#' @return a ggplot
#' @export
#'
#' @examples
#'
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr") %>%
#' pcr_tidy() %>%
#' pcr_plot()
#'


pcr_plot <- function(tidy_pcr) {
        tidy_pcr %>%
                dplyr::filter(!is.na(.data$RQ)) %>%
                ggplot2::ggplot(ggplot2::aes(x = .data$`Sample Name`, y = .data$RQ, fill = .data$`Target Name`)) +
                ggplot2::geom_col() +
                ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$`RQ Min`, ymax = .data$`RQ Max`)) +
                ggplot2::facet_wrap(~.data$`Target Name`, scales = "free") +
                ggplot2::scale_fill_viridis_d(begin = 0.1, end = 0.9) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                      axis.title.x = ggplot2::element_blank())
}
