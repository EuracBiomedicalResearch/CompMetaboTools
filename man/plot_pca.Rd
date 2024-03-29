% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_pca}
\alias{plot_pca}
\title{Plot of PCA results}
\usage{
plot_pca(
  pc,
  pch = 16,
  col = "#000000",
  pc_x = 1,
  pc_y = 2,
  main = "",
  labels = NULL,
  ...
)
}
\arguments{
\item{pc}{the result from a principal component analysis (i.e. the result
returned by \code{prcomp}.}

\item{pch}{the point character. See \code{\link[=plot]{plot()}} or \code{\link[=par]{par()}} for more information.}

\item{col}{the color to be used for each data point/sample.}

\item{pc_x}{\code{integer(1)} defining which principal component should be drawn
on the x-axis.}

\item{pc_y}{\code{integer(1)} defining the principal component to be drawn on the
y-axis.}

\item{main}{\code{character(1)} with the optional title of the plot.}

\item{labels}{\code{character} with length equal to the number of samples. If
provided, these will be displayed instead of data points.}

\item{...}{additional arguments to be passed to the \code{\link[=points]{points()}} or \code{\link[=text]{text()}}
calls (if \code{labels = NULL} or not).}
}
\description{
\code{plot_pca} is a simple utility function to plot the results from a PCA
analysis.
}
\author{
Johannes Rainer
}
