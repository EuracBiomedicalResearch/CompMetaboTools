% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_features.R
\name{groupEicCorrelation}
\alias{groupEicCorrelation}
\title{Group EICs based on their correlation}
\usage{
groupEicCorrelation(
  x,
  aggregationFun = mean,
  threshold = 0.8,
  align = "closest",
  inclusive = FALSE,
  ...
)
}
\arguments{
\item{x}{\code{\link[=MChromatograms]{MChromatograms()}} object of \code{list} of \code{\link[=Chromatogram]{Chromatogram()}} objects.}

\item{aggregationFun}{\code{function} to combine the correlation values between
pairs of EICs across samples (columns). See description for details.}

\item{threshold}{\code{numeric(1)} with the threshold for correlation above which
EICs are grouped together.}

\item{align}{\code{character(1)} defining the method how chromatograms should be
aligned prior correlation. Defaults to \code{align = "closest"}. See
\code{\link[=alignRt]{alignRt()}} for more details.}

\item{inclusive}{\code{logical(1)} defining the grouping approach. With
\code{inclusive = FALSE} (the default) small groups of highly correlated
features are created using the \code{\link[=groupSimilarityMatrix]{groupSimilarityMatrix()}} function. With
\code{inclusive = TRUE} groups are created with features that have at least
one correlation with any other member of the group which is higher than
\code{threshold}.}

\item{...}{parameters for the \code{\link[=correlate]{correlate()}} function for \code{\link[=MChromatograms]{MChromatograms()}}
objects, such as \code{tolerance} to allow specifying the maximal acceptable
difference in retention times between objects. See also \code{\link[=alignRt]{alignRt()}} for
more information.}
}
\value{
\code{factor} same length as \code{nrow(x)} (if \code{x} is a \code{MChromatograms}
object) or \code{length(x)} (if \code{x} is a \code{list}) with the group each EIC
is assigned to.
}
\description{
\code{groupEicCorrelation} groups (extracted ion) chromatograms (EICs) based on
their correlation with each other. If this correlation is \code{>=} than the
provided \code{threshold} they are grouped.

If \code{x} is a \code{\link[=MChromatograms]{MChromatograms()}} object with more than one column (sample),
pairwise correlations between EICs are first calculated for each column
(sample) of \code{x} separately and subsequently aggregated across samples using
\code{aggregationFun}. If \code{x} is a \code{MChromatograms} with 4 rows (EICs) and 3
columns (samples), pairwise correlations are first calculated between all
4 EICs in each of the 3 columns resulting in 3 correlation matrices (of
dimension 4x4). These correlation matrices are combined into a single matrix
by combining the 3 correlation values per comparison with
\code{aggregationFun}. By default the mean of the correlation value between e.g.
EIC 1 and EIC 2 in each of the 3 columns is used as the final correlation
value. Similar to the one-column case EICs are grouped if their (aggregated)
correlation coefficient is larger than \code{threshold}.

Two types of groupings are available:
\itemize{
\item \code{inclusive = FALSE} (the default): the algorithm creates small groups of
highly correlated members, all of which have a correlation with each other
that are \verb{>= threshold}. Note that with this algorithm, rows in \code{x} could
still have a correlation \verb{>= threshold} with one or more elements of a
group they are not part of. See notes below for more information.
\item \code{inclusive = TRUE}: the algorithm creates large groups containing rows that
have a correlation \verb{>= threshold} with at least one element of that group.
For example, if row 1 and 3 have a correlation above the threshold and
rows 3 and 5 too (but correlation between 1 and 5 is below the threshold)
all 3 are grouped into the same group (i.e. rows 1, 3 \strong{and} 5).
}

For more information see \code{\link[=groupByCorrelation]{groupByCorrelation()}}.

Note that it might be useful to set \code{tolerance = 0} if chromatograms from
the \strong{same} sample are compared. This forces retention times of the compared
chromatograms' intensities to be identical.
}
\examples{

library(MSnbase)
set.seed(123)
chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
    intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
    intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
    intensity = c(53, 80, 130, 15, 5, 3, 2))
chrs <- MChromatograms(list(chr1, chr2, chr3))

groupEicCorrelation(chrs)

## With a MChromatograms with two columns, use the maximal correlation
## coefficient found in each of the columns
chrs <- MChromatograms(list(chr1, chr2, chr3, chr1, chr2, chr3), ncol = 2)
groupEicCorrelation(chrs, aggregationFun = max)
}
\seealso{
Other grouping operations: 
\code{\link{groupByCorrelation}()},
\code{\link{groupToSinglePolarityPairs}()}
}
\author{
Johannes Rainer
}
\concept{grouping operations}
