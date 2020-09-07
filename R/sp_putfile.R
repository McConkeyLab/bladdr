#' Upload file to SharePoint
#'
#' @param file path to file to be uploaded
#' @param dest_path path to folder on SharePoint where file is to be uploaded
#' @param token AzureGraph token - see example.
#' @param overwrite should the file overwrite that of one with a similar (see Details) name if it exists?
#'
#' @return an http response
#' @export
#'
#' @details
#' SharePoint API does not care about upper versus lower case.
#' As a matter of good practice, names should be distinct beyond just capitalization.
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
sp_putfile <- function(file, dest_path, token, overwrite = F) {

        base_url <- "https://graph.microsoft.com/v1.0/sites/"

        # Get file extension
        ext_ind <- regexpr("\\..[^\\.]*$", file, perl = T)
        ext     <- regmatches(file, ext_ind)


        # Get file name (with extension)
        file_name_ind <- regexpr("(?<=\\/).[^\\/]*$", file, perl = T)
        file_name     <- regmatches(file, file_name_ind)


        # Parse destination filepath
        parsed <- parse_path(dest_path)


        # Get SharePoint ID
        url    <- paste0(base_url, parsed$host, ":/sites/", parsed$site)
        sp     <- AzureGraph::call_graph_url(token, url, http_verb = "GET")
        sp_id  <- sp$id


        # Get destination folder ID
        url       <- paste0(base_url, sp_id, "/drive/root:/", parsed$rest)
        folder    <- AzureGraph::call_graph_url(token, url, http_verb = "GET")
        folder_id <- folder$id

        # Check if file of similar name exists
        if (!overwrite) {
                url <- paste0(base_url, sp_id, "/drive/items/", folder_id,"/children")
                response <- AzureGraph::call_graph_url(token = token, url = url) %>%
                        dplyr::as_tibble() %>%
                        tidyr::unnest_wider(col = .data$value) %>%
                        dplyr::mutate(name  = tolower(.data$name)) %>%
                        dplyr::filter(.data$name == tolower(file_name))
                if (nrow(response) > 0) stop("\nFile already exists in destination. \nChange filename (case-INsensitive) or set overwrite to TRUE")
        }

        # Upload file
        url <- paste0(base_url, sp_id, "/drive/items/", folder_id,":/", file_name, ":/content")
        AzureGraph::call_graph_url(token     = token,
                                   url       = url,
                                   body      = httr::upload_file(file),
                                   encode    = mime::guess_type(ext),
                                   http_verb = "PUT")
}

