#' Make a heatmap of a select number of genes
#'
#'
#' @param dds a `dds` object with HGCN symbols as rownames
#' @param genes a character vector containing the names of HGCN symbols to be plotted in the heatmap
#' @param stratifier optional. A character vector specifying a column in colData by which to color samples in the pheatmap.
#' @param assay the assay number that contains normalized counts
#' @param ... additional arguments that will be passed to `pheatmap`
#'
#' @return a `pheatmap`
#' @export

eda_heatmap <- function(dds, genes, stratifier = NULL, assay = 2, ...) {
        col_annot <- NULL
        if (!is.null(stratifier)) {
                col_annot <- data.frame(row.names = colnames(dds),
                                        strat = dds[[stratifier]])
        }

        green_black_red <- grDevices::colorRampPalette(c("green", "black", "red"))

        filtered_genes <- dds[which(rownames(dds) %in% genes),]
        norm_exp <- SummarizedExperiment::assay(filtered_genes, assay) %>%
                pheatmap::pheatmap(annotation_col = col_annot, scale = "row", color = green_black_red(10), border_color = NA, ...)
}
