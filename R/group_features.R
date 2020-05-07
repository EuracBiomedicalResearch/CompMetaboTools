#' @description
#'
#' Utility function to group rows/columns in an nxn logical matrix which
#' have a `TRUE`.
#'
#' @param x `logical` `matrix` with same number of rows and columns. See below
#'      for examples.
#'
#' @return `list` of `integer` indexes of rows (columns) that are grouped.
#'
#' @noRd
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Find groups of rows with a correlation higher than 0.9
#' x <- rbind(
#'     c(1, 3, 2, 5),
#'     c(2, 6, 4, 7),
#'     c(1, 1, 3, 1),
#'     c(1, 3, 3, 6),
#'     c(0, 4, 3, 1))
#' xcor <- cor(t(x))
#'
#' .group_logic_matrix(xcor > 0.9)
.group_logic_matrix <- function(x) {
    x <- unname(x)
    nr <- nrow(x)
    if (nr != ncol(x))
        stop("'x' is supposed to be a symmetric matrix")
    groups <- vector("list", nr)
    is_in_group <- rep(FALSE, nr)
    for (i in seq_len(nr)) {
        if (is_in_group[i])
            next
        idx <- which(x[i, ])
        if (length(idx) > 1)
            idx <- unique(c(idx, which(rowSums(x[, idx], na.rm = TRUE) > 0)))
        groups[[i]] <- idx
        is_in_group[idx] <- TRUE
    }
    groups[lengths(groups) > 0]
}

#' @note this expects a `list` of `integer` were each index is only present
#'     **once** and in addition no indices should be missing!
#'
#' @noRd
.index_list_to_factor <- function(x) {
    len <- sum(lengths(x))
    res <- integer(len)
    for (i in seq_along(x))
        res[x[[i]]] <- i
    as.factor(res)
}

#' @title Group rows in a matrix based on their correlation
#'
#' @description
#'
#' The `groupByCorrelation` allows to group rows in a numeric matrix based on
#' their correlation with each other. Note that if for example row 1 and 3 have
#' a correlation above the threshold and rows 3 and 5 too (but correlation
#' between 1 and 5 is below the threshold) all 3 are grouped into the same
#' group (i.e. rows 1, 3 **and** 5).
#'
#' Note that with parameter `f` it is also possible to pre-define groups of
#' rows that should be further sub-grouped based on correlation with each other.
#' In other words, if `f` is provided, correlations are calculated only between
#' rows with the same value in `f` and hence these pre-defined groups of rows
#' are further sub-grouped based on pairwise correlation. The returned `factor`
#' is then `f` with the additional subgroup appended (and separated with a
#' `"."`). See examples below.
#'
#' @param x `numeric` `matrix` where rows should be grouped based on
#'     correlation of their values across columns being larger than `threshold`.
#'
#' @param method `character(1)` with the method to be used for correlation. See
#'     [corr()] for options.
#'
#' @param use `character(1)` defining which values should be used for the
#'     correlation. See [corr()] for details.
#'
#' @param threshold `numeric(1)` defining the cut of value above which
#'     rows are considered to be correlated and hence grouped.
#'
#' @param f optional vector of length equal to `nrow(x)` pre-defining groups
#'     of rows in `x` that should be further sub-grouped. See description for
#'     details.
#'
#' @return `factor` with same length than `nrow(x)` with the group each row
#'     is assigned to.
#'
#' @author Johannes Rainer
#'
#' @importFrom stats cor
#' 
#' @family grouping operations
#'
#' @export
#'
#' @examples
#'
#' x <- rbind(
#'     c(1, 3, 2, 5),
#'     c(2, 6, 4, 7),
#'     c(1, 1, 3, 1),
#'     c(1, 3, 3, 6),
#'     c(0, 4, 3, 1),
#'     c(1, 4, 2, 6),
#'     c(2, 8, 2, 12))
#'
#' ## define which rows have a high correlation with each other
#' groupByCorrelation(x)
#'
#' ## assuming we have some prior grouping of rows, further sub-group them
#' ## based on pairwise correlation.
#' f <- c(1, 2, 2, 1, 1, 2, 2)
#' groupByCorrelation(x, f = f)
groupByCorrelation <- function(x, method = "pearson",
                               use = "pairwise.complete.obs",
                               threshold = 0.9, f = NULL) {
    if (length(threshold) > 1)
        stop("'threshold' has to be of length 1")
    if (!is.null(f)) {
        if (length(f) != nrow(x))
            stop("If 'f' is provided its length has to be equal to 'nrow(x)'")
        f <- as.character(f)
        fnew <- character(length(f))
        for (fg in unique(f)) {
            idx <- which(f == fg)
            idxl <- length(idx)
            if (idxl > 1) {
                cors <- cor(t(x[idx, ]), method = method, use = use) > threshold
                fnew[idx] <- paste0(
                    fg, ".",
                    as.character(.index_list_to_factor(.group_logic_matrix(cors))))
            } else
                fnew[idx] <- paste0(fg, ".1")
        }
        as.factor(fnew)
    } else {
        cors <- cor(t(x), method = method, use = use) > threshold
        .index_list_to_factor(.group_logic_matrix(cors))
    }
}

#' @title Group EICs based on their correlation
#'
#' @description
#'
#' `groupEicCorrelation` groups (extracted ion) chromatograms (EICs) based on
#' their correlation with each other. If this correlation is higher than the
#' provided `threshold` they are grouped.
#'
#' If `x` is a [Chromatograms()] object with more than one column (sample),
#' pairwise correlations between EICs are first calculated for each column
#' (sample) of `x` separately and subsequently aggregated across samples using
#' `aggregationFun`. If `x` is a `Chromatograms` with 4 rows (EICs) and 3
#' columns (samples), pairwise correlations are first calculated between all
#' 4 EICs in each of the 3 columns resulting in 3 correlation matrices (of
#' dimension 4x4). These correlation matrices are combined into a single matrix
#' by combining the 3 correlation values per comparison with
#' `aggregationFun`. By default the mean of the correlation value between e.g.
#' EIC 1 and EIC 2 in each of the 3 columns is used as the final correlation
#' value. Similar to the one-column case EICs are grouped if their (aggregated)
#' correlation coefficient is larger than `threshold`.
#'
#' @param x [Chromatograms()] object of `list` of [Chromatogram()] objects.
#'
#' @param aggregationFun `function` to combine the correlation values between
#'     pairs of EICs across samples (columns). See description for details.
#'
#' @param threshold `numeric(1)` with the threshold for correlation above which
#'     EICs are grouped together.
#'
#' @param ... parameters for the [correlate()] function for [Chromatograms()]
#'     objects.
#'
#' @return `factor` same length as `nrow(x)` (if `x` is a `Chromatograms`
#'     object) or `length(x)` (if `x` is a `list`) with the group each EIC
#'     is assigned to.
#' 
#' @importMethodsFrom xcms correlate
#'
#' @importFrom MSnbase Chromatograms
#' 
#' @family grouping operations
#'
#' @author Johannes Rainer
#' 
#' @export
#'
#' @examples
#'
#' library(MSnbase)
#' set.seed(123)
#' chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
#'     intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
#' chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
#'     intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
#' chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
#'     intensity = c(53, 80, 130, 15, 5, 3, 2))
#' chrs <- Chromatograms(list(chr1, chr2, chr3))
#'
#' groupEicCorrelation(chrs)
#'
#' ## With a Chromatograms with two columns, use the maximal correlation
#' ## coefficient found in each of the columns
#' chrs <- Chromatograms(list(chr1, chr2, chr3, chr1, chr2, chr3), ncol = 2)
#' groupEicCorrelation(chrs, aggregationFun = max)
groupEicCorrelation <- function(x, aggregationFun = mean,
                                threshold = 0.8, ...) {
    nr <- nrow(x)
    nc <- ncol(x)
    res <- array(NA_real_, dim = c(nr, nr, nc))
    ## For performance issues it would also be possible to run with full = FALSE
    for (i in seq_len(nc)) {
        res[, , i] <- correlate(x[, i], ...)
    }
    res <- apply(res, c(1, 2), aggregationFun, na.rm = TRUE) > threshold
    .index_list_to_factor(.group_logic_matrix(res))
}
