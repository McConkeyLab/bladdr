#' Summon cell line data
#' @details This function requires access to the GBCI SharePoint.
#' @return a DDS containing cell-line RNA expression
#' @export
summon_cell_data <- function() {
  get_gbci_file("Datasets/Cell Line RNA Seq Data/thirty-cell-lines-rnaseq/cell-lines_norm_clades.Rds") |>
    readRDS()
}
