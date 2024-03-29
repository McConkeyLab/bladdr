% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dilute}
\alias{dilute}
\title{Calculate dilution from known concentrations}
\usage{
dilute(c1, c2 = min(c1), v2, round_for_pipettes = TRUE, quiet = FALSE)
}
\arguments{
\item{c1}{Numeric. Initial concentration of sample.}

\item{c2}{Numeric. Target concentration of sample.}

\item{v2}{Numeric. Target final volume of sample. If \code{round_for_pipettes = TRUE}, assumes volume is mL.}

\item{round_for_pipettes}{Logical. If TRUE, rounds values to the accuracy of
standard pipettes using \code{make_pipette_vol}.}

\item{quiet}{Logical. If FALSE, will warn when dilution is impossible to do
without concentrating sample.}
}
\value{
a named list, with \code{sample_to_add} as the volume of sample to add,
and \code{add_to} as the volume to dilute the sample into.
}
\description{
Calculate dilution from known concentrations
}
\examples{
dilute(203, 70, 10)
dilute(203, 70, 10, round_for_pipettes = FALSE)

}
