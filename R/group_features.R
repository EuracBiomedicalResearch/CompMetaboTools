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
#' @return `list` of indices of rows which are grouped based on their
#'     correlation.
#'
#' @author Johannes Rainer
#'
#' @importFrom stats cor
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
#'     c(0, 4, 3, 1))
#'
#' groupByCorrelation(x)
groupByCorrelation <- function(x, method = "pearson",
                               use = "pairwise.complete.obs",
                               threshold = 0.9) {
    if (length(threshold) > 1)
        stop("'threshold' has to be of length 1")
    cors <- cor(t(x), method = method, use = use) > threshold
    .group_logic_matrix(cors)
}
