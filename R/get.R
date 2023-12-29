# A series of thin wrappers around get_gbci for particular types of data
get_f <- function(sub_dir) {
  function(name, user_name = NULL, ...) {
    if (is.null(user_name)) user_name <- get_username_from_env()
    path <- fs::path("Raw Data", sub_dir, user_name, name)
    get_gbci(path = path, ...)
  }
}

#' Get from a raw data folder
#'
#' A collection of thin wrappers around `get_gbci`
#'
#' @param name Character, the name of the file
#' @param user_name Optional character. Should be the name used for your
#'   personal data directories in the SharePoint. If NULL, will use the
#'   GBCI_USERNAME environmental variable.
#' @param ... Additional arguments passed to `get_gbci`
#'
#' @export
get_pcr <- get_f("qPCR")

#' @export
get_nanodrop <- get_f("Nanodrop")

#' @export
get_chemidoc <- get_f("ChemiDoc")

#' @export
get_licor <- get_f("licor")

#' @export
get_evos <- get_f("EVOS")

#' @export
get_qubit <- get_f("Qubit")

#' @export
get_spectramax <- get_f("SPECTRAmax")
