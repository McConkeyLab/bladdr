
#' @noRd
get_gene_index <- function(dds, gene) {
        ind <- which(rownames(dds) == gene)
        if(identical(ind, integer(0))) {
                stop("Gene does not exist in this dataset")
        }
        ind
}

theme_tufte <- function(font_size = 30) {

        if (Sys.info()[["sysname"]] == "Windows") {
                font <- "Gill Sans MT"
        } else if (Sys.info()[["sysname"]] == "Darwin") {
                font <- "Gill Sans"
        }

        ggplot2::theme(
                panel.grid = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = "#FFFFF8", color = "#CCCCCC"),
                plot.background = ggplot2::element_rect(fill = "#FFFFF8"),
                strip.background = ggplot2::element_rect(fill = "#BBBBB0"),
                legend.background = ggplot2::element_rect(fill = "#FFFFF8"),
                legend.position = "top",
                legend.key = ggplot2::element_blank(),
                text = ggplot2::element_text(family = font, size = font_size)
        )
}
