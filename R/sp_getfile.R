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
#' gr <- AzureGraph::create_graph_login() # Will trigger a login
#' me <- gr$get_user("me")
#' token <- me$token}
#'
sp_getfile <- function(path, token) {

        tmp    <- tempfile()
        parsed <- parse_path(path)
        url    <- paste0("https://graph.microsoft.com/v1.0/sites/", parsed$host, ":/sites/", parsed$site)
        sp_dat <- AzureGraph::call_graph_url(token, url, http_verb = "GET")
        sp_id  <- sp_dat$id
        file_path <- regexpr(paste0("(?<=", parsed$drive, "/).*$"), path, perl = T)
        file_path <- regmatches(path, file_path)
        url <- paste0("https://graph.microsoft.com/v1.0/sites/", sp_id, "/drive/root:/", file_path, ":/content")
        file <- AzureGraph::call_graph_url(token, url,
                                 httr::write_disk(tmp, overwrite = T),
                                 http_verb = "GET")
        tmp
}
