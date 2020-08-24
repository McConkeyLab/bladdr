pcr_plate_view <- function(tidy_pcr) {
        tidy_pcr <- tidy_pcr %>%
                dplyr::mutate(well_row = stringr::str_extract(`Well Position`, "^.{1}"),
                              well_col = as.numeric(stringr::str_extract(`Well Position`, "[:digit:]{1,2}$")),
                              well_row = as.numeric(factor(well_row, levels = LETTERS)))

        ggplot2::ggplot(tidy_pcr, ggplot2::aes(x = well_col, y = well_row, fill = `Target Name`)) +
                ggplot2::geom_tile(ggplot2::aes(size = 2)) +
                ggplot2::coord_cartesian(xlim = c(1,24), ylim = c(16, 1)) +
                ggplot2::scale_y_continuous(breaks = 1:16, labels = LETTERS[1:16]) +
                ggplot2::scale_x_continuous(breaks = 1:24, minor_breaks = NULL) +
                ggplot2::guides(size = F)
}
