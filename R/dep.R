#' Quantify protein concentration from a MicroBCA assay
#'
#' @param x A `spectramax`, `gp`, or `data.frame` object, or path to SPECTRAmax
#' .xls(x)/.txt file.
#' @param replicate_orientation Either 'h' or 'v' - see Details.
#' @param sample_names Character vector of sample names.
#' @param remove_empty Should wells that have less absorbance than the lowest
#' standard be dropped?
#'
#' @details The standards must be in ascending concentration starting in the
#'   upper left corner. Whether this is from from left to right or top to bottom
#'   can be specified in 'replicate orientation'. Note that 'replicate
#'   orientation' specified the direction that REPLICATES lie, NOT the direction
#'   the samples flow (which will be opposite).
#'
#' @return a `tibble`
#' @export
qp <- function(x,
               replicate_orientation = c("h", "v"),
               sample_names = NULL,
               remove_empty = TRUE,
               remove_outliers = c("all", "samples", "standards", "none")) {
  lifecycle::deprecate_stop("1.0.0", "qp()", "qp::qp()")
}

# Fit conc ~ abs using standards absorbances -----------------------------------
qp_fit <- function(x) {
  lifecycle::deprecate_stop("1.0.0", "qp_fit()", "qp::qp_fit()")
}

# Calculate Dilutions ----------------------------------------------------------
#' Calculate dilutions for an analyzed `qp` `list`
#'
#' @param x The output of `qp()`
#' @param target_conc Target concentration in (mg/mL) protein
#' @param target_vol Target volume in uL
#'
#' @return A list, where the `qp` item has volumes of lysate and volumes of H2O
#' to add.
#' @export
qp_calc_dil <- function(x, target_conc, target_vol) {
  lifecycle::deprecate_stop("1.0.0", "qp_calc_dil()", "qp::qp_dilute()")
}

# Visualization ----------------------------------------------------------------
#' View the absorbances of an analyzed `qp` as they were on the plate
#'
#' @param x The output of `qp()`
#'
#' @return a `ggplot`
#' @export
make_qp_plate_view <- function(x) {
  lifecycle::deprecate_stop(
    "1.0.0", "make_qp_plate_view()", "qp::qp_plot_plate()"
  )
}

#' View an absorbance/concentration plot
#'
#' @param x The output of `qp()`
#'
#' @return a `ggplot`
#' @export
make_qp_standard_plot <- function(x) {
  lifecycle::deprecate_stop(
    "1.0.0", "make_qp_standard_plot()", "qp::qp_plot_standards()"
  )
}

#' Read in data into a common quantify protein format
#' @param x A `gp`, `data.frame`/`tibble`, or character
#' path to a raw SPECTRAmax .xls(x)/.txt
#' @param ... Unused
#'
#' @return A `gp`
#' @export
qp_read <- function(x, ...) {
  lifecycle::deprecate_stop("1.0.0", "qp_read()", "qp::qp_tidy()")
}
