% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MChromatograms.R
\name{c,MChromatograms-method}
\alias{c,MChromatograms-method}
\alias{c,XChromatograms-method}
\title{Combine MChromatograms objects}
\usage{
\S4method{c}{MChromatograms}(x, ...)

\S4method{c}{XChromatograms}(x, ...)
}
\arguments{
\item{x}{\code{MChromatograms} object.}

\item{...}{\code{MChromatograms} object that should be appended to \code{x}.}
}
\value{
\code{MChromatograms} object.
}
\description{
\code{\link[=MChromatograms]{MChromatograms()}} objects can be concatenated (row-wise) with the \code{c}
function resulting in a \code{MChromatograms} object with the same number of
samples (columns).
}
\note{
The \code{c} function for \code{XChromatograms} objects (defined in the \code{xcms} package)
is not supported.
}
\author{
Johannes Rainer
}
