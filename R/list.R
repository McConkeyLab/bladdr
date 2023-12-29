list_f <- function(sub_dir) {
  function(user_name = NULL, ...) {
    if (is.null(user_name)) user_name <- get_username_from_env()
    path <- fs::path("Raw Data", sub_dir, user_name)
    list_gbci_dir(path = path, ...)
  }
}

#' List from a raw data folder
#'
#' A collection of thin wrappers around `list_gbci_dir`.
#'
#' @param user_name Optional character. Should be the name used for your
#'   personal data directories in the SharePoint. If NULL, will use the
#'   GBCI_USERNAME environmental variable.
#' @param ... Additional arguments passed to `list_gbci_dir`
#'
#' @export
list_pcr <- list_f("qPCR")

#' @rdname list_pcr
#' @export
list_nanodrop <- list_f("Nanodrop")

#' @rdname list_pcr
#' @export
list_chemidoc <- list_f("ChemiDoc")

#' @rdname list_pcr
#' @export
list_licor <- list_f("licor")

#' @rdname list_pcr
#' @export
list_evos <- list_f("EVOS")

#' @rdname list_pcr
#' @export
list_qubit <- list_f("Qubit")

#' @rdname list_pcr
#' @export
list_spectramax <- list_f("SPECTRAmax")
