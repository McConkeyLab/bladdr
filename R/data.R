#' Example data of RNA samples with concentrations
#'
#' A dataset containing fabricated sample names and RNA concentrations
#'
#' @format A data frame with 8 rows and 2 columns
#' \describe{
#'   \item{sample}{name of sample}
#'   \item{conc}{concentration of RNA, in ng/uL}
#' }
"dummy_rna_conc"

#' Example data from an MTT
#'
#' A dataset including absorbance data at 562nm and 660nm
#'
#' There are four conditions on this plate. The plate is divided into quadrants.
#' Concentrations of drug increase from left to right (0, 1nM, 10nM, 100nM, 1uM,
#' 10uM).
#'
#' @format A `spectramax` object
#' \describe{
#'   \item{data}{a `gp` object containing absorbance values}
#' }
"mtt"
