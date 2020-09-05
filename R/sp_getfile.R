#' Download File From SharePoint
#'
#' @param path a character vector. Obtained by going online to SharePoint, right clicking on file, clicking 'details', clicking 'more details' in the resultant side pane, then clicking the 'copy' button right next to the "Path" header.
#' @param token an AzureGraph token. See Examples.
#'
#' @return a tempfile path pointing to the location of the file downloaded
#' @export
#'
#' @examples
#'
#' # Getting token:
#' \dontrun{
#' gr <- create_graph_login() # Will trigger a login
#' me <- gr$get_user("me")
#' token <- me$token
#' }
#'
sp_getfile <- function(path, token) {

        tmp <- tempfile()

        path_trimmed <- stringr::str_remove(path, "https://livejohnshopkins.sharepoint.com/sites/")
        sp <- stringr::str_extract(path_trimmed, "^[^/]*")
        url <- paste0("https://graph.microsoft.com/v1.0/sites/livejohnshopkins.sharepoint.com:/sites/", sp)
        sp_dat <- AzureGraph::call_graph_url(token, url, http_verb = "GET")
        sp_id <- sp_dat$id

        path_dbl_trimmed <- stringr::str_remove(path_trimmed, paste0(sp, "/Shared%20Documents/"))
        url <- paste0("https://graph.microsoft.com/v1.0/sites/", sp_id, "/drive/root:/", path_dbl_trimmed, ":/content")
        file <- AzureGraph::call_graph_url(token, url,
                                httr::write_disk(tmp, overwrite = T),
                                http_verb = "GET")
        tmp
}
