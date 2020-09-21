#' View sample plating layout
#'
#' @param tidy_pcr an output from the `pcr_tidy` function
#' @param fill the variable to use to fill the geom_tiles
#'
#' @return a ggplot
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#'
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr") %>%
#' pcr_tidy() %>%
#' pcr_plate_view()
#'
pcr_plate_view <- function(tidy_pcr, fill = .data$target_name) {
        usr_fill <- substitute(fill)
        tidy_pcr %>%
                if (tidy_pcr$plate_type == "384-Well Block"){
                ggplot2::ggplot(ggplot2::aes(x = .data$well_col, y = .data$well_row, fill = eval(usr_fill))) +
                ggplot2::geom_tile(ggplot2::aes(size = 2)) +
                ggplot2::coord_cartesian(xlim = c(1,24), ylim = c(16, 1)) +
                ggplot2::scale_y_continuous(breaks = 1:16, labels = LETTERS[1:16]) +
                ggplot2::scale_x_continuous(breaks = 1:24, minor_breaks = NULL) +
                ggplot2::labs(fill = deparse(usr_fill)) +
                ggplot2::guides(size = F)}
        else if (tidy_pcr$plate_type == "96-Well Block") {
                ggplot2::ggplot(ggplot2::aes(x = .data$well_col, y = .data$well_row, fill = eval(usr_fill))) +
                        ggplot2::geom_tile(ggplot2::aes(size = 2)) +
                        ggplot2::coord_cartesian(xlim = c(1,12), ylim = c(8, 1)) +
                        ggplot2::scale_y_continuous(breaks = 1:8, labels = LETTERS[1:8]) +
                        ggplot2::scale_x_continuous(breaks = 1:12, minor_breaks = NULL) +
                        ggplot2::labs(fill = deparse(usr_fill)) +
                        ggplot2::guides(size = F)}
}
