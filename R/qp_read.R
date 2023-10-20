#' Read in data into a common quantify protein format
#' @param x A `gp`, `data.frame`/`tibble`, or character
#' path to a raw SPECTRAmax .xls(x)/.txt
#' @param ... Unused
#'
#' @return A `gp`
#' @export
qp_read <- function(x, ...) {
  UseMethod("qp_read")
}

#' @export
#' @rdname qp_read
qp_read.character <- function(x, ...) {
  mop::read_spectramax(x, ...) |> qp_read.spectramax()
}

#' @export
#' @rdname qp_read
qp_read.data.frame <- function(x, ...) {
  gplate::as_gp(x)
}

#' @export
#' @rdname qp_read
qp_read.gp <- function(x, ...) {
  x
}

#' @export
#' @rdname qp_read
qp_read.spectramax <- function(x, ...) {
  if (!(562 %in% x$wavelengths)) {
    rlang::warn("x$wavelengths does not contain 562")
  }

  plate_index <- which(sapply(x$data, \(x) x$type == "Plate"))

  if (length(plate_index) != 1) {
    rlang::abort("Supplied data does not include 1 (and only 1) plate")
  }

  gp <- x$data[[plate_index]]$data

  if ("value" %in% colnames(gplate::well_data(gp))) {
    rlang::abort("SPECTRAmax has an unexpected 'value' column")
  } else {
    wd <- gplate::well_data(gp)
    wd$value <- wd$nm562
    gp$well_data <- wd
  }

  gp
}
