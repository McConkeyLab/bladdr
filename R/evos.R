#' Flatten an EVOS timelapse directory
#'
#' @param path character. Path to directory that contains all timelapses
#'
#' @return Invisibly returns logical vector the length of the number of files.
#'   As a side effect, all files in that directory will be copied to the top
#'   level of that directory with a file name that is a concatination of all
#'   parent directories, separated by an underscore.
#' @export
evos_flatten <- function(path, parse_names = FALSE) {
  path <- gsub("/$", "", path)
  file_paths <- list.files(path, recursive = TRUE, full.names = TRUE, pattern = "*.tif")
  new_file_paths <- sapply(strsplit(file_paths, "/"), \(x)
                           paste0(
                             paste(head(x, -3), collapse = "/"),
                             "/",
                             paste(tail(x, 3), collapse = "_")))
  if (parse_names) new_file_paths <- lapply(new_file_paths, evos_parse_name) |> unlist()
  file.copy(file_paths, new_file_paths, overwrite = FALSE)
}


#' Sensibly parse default EVOS file names
#'
#' @param names
#'
#' @return A filename that takes the form of
#'   "parent-directory/b{N}_{YYYY-MM-DD}<_filter>.tif" Where {N} is the number
#'   of the beacon and <_filter> is an underscore prepended filter name (eg
#'   _phase) which may or may not exist.
#' @noRd
evos_parse_name <- function(name) {
  fname <-  basename(name)
  fname_split <- strsplit(fname, split = " ")[[1]]
  date <- fname_split[4]
  flattened_split <- strsplit(fname_split[5], split = "_")[[1]]
  beacon <- paste0("b", strsplit(flattened_split[2], "-")[[1]][2])
  if (length(flattened_split) > 3)
    return(paste0(dirname(name), "/", beacon, "_", date, "_", tolower(flattened_split[4])))
  paste0(dirname(name), "/", beacon, "_", date, ".tif")
}
