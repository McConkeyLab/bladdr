#' Plot the expression of one gene across a stratifying argument
#'
#' @param dds A `DESeqDataSet` object
#' @param gene A character vector HGNC Symbol, like "KRT14"
#' @param stratifier A name of a column of the `dds` `colData` that will be used as the x-axis to stratify gene expression into groups
#' @param assay The index of the assay that contains the normalized counts. Defaults to `2`.
#'
#' @return A ggplot
#' @export
#'
#'
eda_one_gene <- function(dds, gene, stratifier, assay = 2) {

        gene_ind <- which(rownames(dds) == gene)

        if(identical(gene_ind, integer(0))) {
                stop("Gene does not exist in this dataset")
        }

        expression <- SummarizedExperiment::assay(dds, assay)[gene_ind,]
        plotting_data = data.frame(as.data.frame(SummarizedExperiment::colData(dds)), expression = expression)
        ggplot2::ggplot(plotting_data, ggplot2::aes(x = {{stratifier}}, y = expression, color = {{stratifier}}, fill = {{stratifier}})) +
                ggplot2::geom_jitter(alpha = 0.7) +
                ggplot2::labs(y = gene, x = deparse(substitute(stratifier))) +
                ggplot2::theme_minimal() +
                ggplot2::scale_color_viridis_d(option = "B", begin = 0.5, end = 0.8) +
                ggplot2::scale_fill_viridis_d(option = "B", begin = 0.5, end = 0.8) +
                ggplot2::coord_cartesian(ylim = c(0, NA)) +
                ggplot2::theme(panel.grid = ggplot2::element_blank(),
                      text = ggplot2::element_text(size = 20),
                      legend.position = "none",
                      panel.background = ggplot2::element_rect(fill = "#FFFFF8", color = NA))
}
