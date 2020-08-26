#' View sample plating layout
#'
#' @param tidy_pcr an output from the `pcr_tidy` function
#'
#' @return a ggplot
#' @export
#'
#' @examples
#'
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr") %>%
#' pcr_tidy() %>%
#' pcr_plate_view()
#'
pcr_plate_view <- function(tidy_pcr) {
        tidy_pcr %>%
                ggplot2::ggplot(ggplot2::aes(x = .data$well_col, y = .data$well_row, fill = .data$`Target Name`)) +
                ggplot2::geom_tile(ggplot2::aes(size = 2)) +
                ggplot2::coord_cartesian(xlim = c(1,24), ylim = c(16, 1)) +
                ggplot2::scale_y_continuous(breaks = 1:16, labels = LETTERS[1:16]) +
                ggplot2::scale_x_continuous(breaks = 1:24, minor_breaks = NULL) +
                ggplot2::guides(size = F)
}
