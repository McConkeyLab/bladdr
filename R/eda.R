eda_two_gene <- function(dds, gene1, gene2, color_by_clade = TRUE) {
        get_gene_ind <- function(gene_name) {

        }
}


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

        gene_ind <- get_gene_index(dds, gene)
        expression <- SummarizedExperiment::assay(dds, assay)[gene_ind,]
        plotting_data = data.frame(as.data.frame(SummarizedExperiment::colData(dds)), expression = expression)
        ggplot2::ggplot(plotting_data, ggplot2::aes(x = .data[[stratifier]], y = expression, color = .data[[stratifier]], fill = .data[[stratifier]])) +
                ggplot2::geom_dotplot(binaxis = "y", binwidth = 0.2, stackdir = "center") +
                ggplot2::labs(y = gene, x = NULL) +
                ggrepel::geom_text_repel(aes(label = .data$cell), nudge_x = 0.5) +

}
