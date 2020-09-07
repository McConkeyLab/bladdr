#' Upload file to SharePoint
#'
#' @param file path to file to be uploaded
#' @param dest_path path to folder on SharePoint where file is to be uploaded
#' @param token AzureGraph token - see example.
#'
#' @return an http response
#' @export
#'
#' @examples
#'
#' # Getting token:
#' \dontrun{
#'
#' gr <- AzureGraph::create_graph_login() # Will trigger a login
#' me <- gr$get_user("me")
#' token <- me$token
#' }
sp_putfile <- function(file, dest_path, token) {


        # Get file extension
        ext_ind <- regexpr("\\..[^\\.]*$", file, perl = T)
        ext     <- regmatches(file, ext_ind)


        # Get file name (with extension)
        file_name_ind <- regexpr("(?<=\\/).[^\\/]*$", file, perl = T)
        file_name     <- regmatches(file, file_name_ind)


        # Parse destination filepath
        parsed <- parse_path(dest_path)


        # Get SharePoint ID
        url    <- paste0("https://graph.microsoft.com/v1.0/sites/", parsed$host, ":/sites/", parsed$site)
        sp     <- AzureGraph::call_graph_url(token, url, http_verb = "GET")
        sp_id  <- sp$id


        # Get Destination Folder ID
        url       <- paste0("https://graph.microsoft.com/v1.0/sites/", sp_id, "/drive/root:/", parsed$rest)
        folder    <- AzureGraph::call_graph_url(token, url, http_verb = "GET")
        folder_id <- folder$id


        # Upload File
        url <- paste0("https://graph.microsoft.com/v1.0/sites/", sp_id, "/drive/items/", folder_id,":/", file_name, ":/content")
        AzureGraph::call_graph_url(token = token,
                                   url = url,
                                   body = httr::upload_file(file),
                                   encode = mime::guess_type(ext),
                                   http_verb = "PUT")

}
