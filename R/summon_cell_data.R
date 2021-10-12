#' Summon cell line data
#' @details This function requires access to the GBCI SharePoint.
#' @return a DDS containing cell-line RNA expression
#' @export
summon_cell_data <- function() {
        cell_line_data_url <- "https://livejohnshopkins.sharepoint.com/sites/GBCIStorage/Shared%20Documents/Datasets/Cell%20Line%20RNA%20Seq%20Data/thirty-cell-lines-rnaseq/cell-lines_norm_clades.Rds"
        token <- sharepointr::sharepoint_token()
        cell_df <- sharepointr::sharepoint_get(cell_line_data_url, token = token) |>
                readRDS()
}

