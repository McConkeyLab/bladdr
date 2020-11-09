#' View sample plating layout
#'
#' @param tidy_pcr an output from the `pcr_tidy` function
#' @param fill the variable to use to fill the geom_tiles
#'
#' @return a ggplot
#' @export
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 aes
#'
#' @examples
#'
#' dat_path <- system.file("extdata", "untidy-pcr-example.xls", package = "bladdr") %>%
#' pcr_tidy() %>%
#' pcr_plate_view()

pcr_plate_view <- function(tidy_pcr, fill = .data$target_name) {
        usr_fill <- substitute(fill)
        if (sum(grepl("384", tidy_pcr$plate_type), na.rm = T) > 0) {
                x <- 2
        } else if (sum(grepl("96", tidy_pcr$plate_type), na.rm = T) > 0) {
                x <- 1
        }
        ggplot2::ggplot(tidy_pcr, aes(x = .data$well_col, y = .data$well_row, fill = eval(usr_fill))) +
                ggplot2::geom_tile(aes(size = 2)) +
                ggplot2::coord_cartesian(xlim = c(1,(x*12)), ylim = c((x*8), 1)) +
                ggplot2::scale_y_continuous(breaks = 1:(x*8), labels = LETTERS[1:(x*8)]) +
                ggplot2::scale_x_continuous(breaks = 1:(x*12), minor_breaks = NULL) +
                ggplot2::labs(fill = deparse(usr_fill)) +
                ggplot2::guides(size = F)
}
