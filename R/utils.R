#' Get indices of genes by name in a dds
#' @noRd
get_gene_index <- function(dds, genes) {
  inds <- which(rownames(dds) %in% genes)
  if(identical(inds, integer(0))) {
    stop("None of the supplied genes exist in this dataset.")
  } else if (length(genes) != length(inds)) {
    found <- rownames(dds)[inds]
    missing <- setdiff(genes, found)
    warning("Dataset does not contain gene: ", paste(missing, collapse = ", "), "\n")
  }
  inds
}

#' Get expression of gene by indices in a dds
#'
#' @param dds A DESeqDataSet
#' @param gene_indices integer vector that refer to the row number of genes to retrieve
#' @param assay integer, assay index that refers to the desired expression data slot
#'
#' @return A dataframe where row names are sample names, and column names are
#'   gene names
#'
#' @details Importantly, returns data as a `data.frame` without dropping
#'   dimensions.
get_gene_expression <- function(dds, gene_indices, assay) {
  SummarizedExperiment::assay(dds, assay)[gene_indices, , drop = FALSE] |>
    t()
}


#' ggplot theme in the style of Edward Tufte
#'
#' @param font_size numeric. Size of `element_text` font.
#'
#' @return A ggplot theme
#' @noRd
theme_tufte <- function(font_size = 30) {

  if (Sys.info()[["sysname"]] == "Windows") {
    font <- "Gill Sans MT"
  } else if (Sys.info()[["sysname"]] == "Darwin") {
    font <- "Gill Sans"
  } else {
    font <- "Arial"
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
