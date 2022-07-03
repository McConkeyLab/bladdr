#' Summon cell line data
#' @details This function requires access to the GBCI SharePoint.
#' @return a DDS containing cell-line RNA expression
#' @export
summon_cell_data <- function() {
  sp <- Microsoft365R::get_sharepoint_site(site_url = "https://livejohnshopkins.sharepoint.com/sites/GBCIStorage")
  drive <- sp$get_drive()
  temp <- tempfile()
  drive$download_file("Datasets/Cell Line RNA Seq Data/thirty-cell-lines-rnaseq/cell-lines_norm_clades.Rds", dest = temp)
  readRDS(temp)
}

