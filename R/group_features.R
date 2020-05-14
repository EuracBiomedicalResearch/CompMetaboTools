#' @description
#'
#' Utility function to group rows/columns in an nxn logical matrix which
#' have a `TRUE`. Note that this groups also rows/columns together that are not
#' direclty *linked* with a `TRUE`, but that are linked via other rows/columns
#' they have in common (i.e. if between row 2 and 4 is a `TRUE`, but also
#' between 3 and 4, all 3 of them are joined together, even if they are not
#' directly linked).
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
#'
#' xcor <- matrix(FALSE, ncol = 13, nrow = 13)
#' for (i in 1:13)
#'     xcor[i, i] <- TRUE
#' xcor[8, 6] <- TRUE
#' xcor[8, 7] <- TRUE
#' xcor[9, 7] <- TRUE
#' xcor[11, 7] <- TRUE
#' xcor[6, 8] <- TRUE
#' xcor[7, 8] <- TRUE
#' xcor[10, 8] <- TRUE
#' xcor[13, 8] <- TRUE
#' xcor[7, 9] <- TRUE
#' xcor[8, 10] <- TRUE
#' xcor[7, 11] <- TRUE
#' xcor[12, 11] <- TRUE
#' xcor[11, 12] <- TRUE
#' xcor[8, 13] <- TRUE
#' .group_logic_matrix(xcor)
#'
#' xcor <- matrix(FALSE, ncol = 10, nrow = 10)
#' for (i in seq_len(ncol(xcor))) {
#'     xcor[i, i] <- TRUE
#' }
#' xcor[1, 4] <- TRUE
#' xcor[4, 1] <- TRUE
#' xcor[2, 8] <- TRUE
#' xcor[8, 2] <- TRUE
#' xcor[3, 9] <- TRUE
#' xcor[9, 3] <- TRUE
#' xcor[8, 9] <- TRUE
#' xcor[9, 8] <- TRUE
#' .group_logic_matrix(xcor)
.group_logic_matrix <- function(x) {
    x <- unname(x)
    nr <- nrow(x)
    if (nr != ncol(x))
        stop("'x' is supposed to be a symmetric matrix")
    groups <- vector("list", nr)
    idx <- which(x, arr.ind = TRUE)
    is_in_group <- rep(FALSE, nr)
    for (i in seq_len(nr)) {
        if (is_in_group[i])
            next
        elms <- idx[idx[, 1] == i, 2]
        if (length(elms) == 1) {
            groups[[i]] <- i
            is_in_group[i] <- TRUE
        } else {
            while(!all(is_in_group[elms])) {
                is_in_group[elms] <- TRUE
                ## Get all rows containing these elements
                elms <- unique(idx[idx[, 1] %in% elms, 2])
            }
            groups[[i]] <- elms
        }
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
        if (!is.factor(f))
            f <- factor(f, levels = unique(f))
        fnew <- rep(NA_character_, length(f))
        for (fg in levels(f)) {
            idx <- which(f == fg)
            idxl <- length(idx)
            if (idxl > 1) {
                cors <- cor(t(x[idx, ]), method = method, use = use) > threshold
                ## Ensure diagonal matrix is TRUE so that even if some features
                ## have a correlation value of NA they are not dropped
                cors[cbind(1:idxl, 1:idxl)] <- TRUE
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
    ## Ensure diagonal is always TRUE to not drop any features!
    res[cbind(1:nr, 1:nr)] <- TRUE
    .index_list_to_factor(.group_logic_matrix(res))
}

#' @title Sub-group allowing only single positive/polarity pairs per group
#'
#' @description
#'
#' Based on a grouping `f` and the polarity of each element `polarity`, the
#' function ensures each feature group to consist of only a single
#' positive/negative feature pair. Thus, each of the groups in `f` is further
#' sub-grouped into positive/negative feature pairs with the highest correlation
#' of feature values across samples. In the returned grouping no group will
#' have more than one feature of the same polarity.
#'
#' This function can be helpful for merging feature grouping results from
#' positive and negative polarity.
#'
#' @param f vector defining the initial grouping of elements (e.g. such as
#'     returned by `groupByCorrelation`.
#'
#' @param polarity vector (same length than `f`) with the polarity for each
#'     element (feature).
#'
#' @param fvals `numeric` `matrix` (number of rows matching `length(f)`) with
#'     feature values across e.g. samples. Correlations will be calculated
#'     between rows of this matrix.
#'
#' @return `factor` with the grouping, ensuring that each group consists of only
#'     a single positive/negative pair.
#'
#' @author Johannes Rainer
#'
#' @family grouping operations
#'
#' @export
#'
#' @examples
#'
#' ## Define a simple matrix where the first 4 and the last 3 rows are highly
#' ## correlated with each other.
#' x <- rbind(
#'     c(4, 3, 5, 1),
#'     c(4, 2, 5, 1),
#'     c(4, 3, 4, 1),
#'     c(4, 3, 4, 1),
#'     c(4, 4, 4, 9),
#'     c(4, 4, 4, 9),
#'     c(4, 4, 4, 9))
#'
#' ## Determin the expected grouping by correlation
#' grp <- groupByCorrelation(x)
#' grp
#'
#' ## Each second row represents however a feature measured in a different
#' ## polarity. From another grouping (e.g. based on EIC correlation) we know
#' ## however that none of the positive polarity features should be
#' ## grouped together. We apply now the function to further subgroup `grp`
#' ## keeping only single positive/negative pairs in each group of `grp`.
#' pol <- c("NEG", "POS", "NEG", "POS", "NEG", "POS", "NEG")
#'
#' groupToSinglePolarityPairs(grp, pol, x)
groupToSinglePolarityPairs <- function(f, polarity = rep(1, length(f)), fvals) {
    if (missing(f) || missing(polarity) || missing(fvals))
        stop("Parameters 'f', 'polarity' and 'fvals' are required")
    if (length(f) != length(polarity))
        stop("Length of 'f' and 'polarity' has to match")
    if (length(f) != nrow(fvals))
        stop("Number of rows of 'fvals' and length of 'f' has to match")
    if (!is.factor(f))
        f <- factor(f, levels = unique(f))
    res <- rep(NA_character_, length(f))
    for (fg in levels(f)) {
        idx <- which(f == fg)
        idxl <- length(idx)
        if (idxl > 1) {
            pols <- polarity[idx]
            if (length(unique(pols)) > 1) {
                cors <- cor(t(fvals[idx, ]), use = "pairwise.complete.obs",
                            method = "pearson")
                cors[!upper.tri(cors)] <- NA
                grpl <- list()
                while(any(!is.na(cors))) {
                    idx_max <- which(cors == max(cors, na.rm = TRUE),
                                     arr.ind = TRUE)[1, , drop = FALSE]
                    cors[idx_max] <- NA
                    idx_max <- as.numeric(idx_max)
                    ## Only consider pairing if we have pos AND neg polarity
                    if (length(unique(pols[idx_max])) > 1) {
                        grpl[[(length(grpl) + 1)]] <- idx_max
                        cors[idx_max, ] <- NA
                        cors[, idx_max] <- NA
                    }
                }
                grps_tmp <- rep(NA_integer_, idxl)
                for (i in seq_along(grpl))
                    grps_tmp[grpl[[i]]] <- i
                nas <- is.na(grps_tmp)
                if (any(nas))
                    grps_tmp[nas] <- seq((length(grpl) + 1),
                                         length.out = sum(nas))
                res[idx] <- paste0(fg, ".", grps_tmp)
            } else
                res[idx] <- paste0(fg, ".", seq_len(idxl))
        } else
            res[idx] <- paste0(fg, ".1")
    }
    as.factor(res)
}

#' @title Grouping of values into sets with smallest differences
#'
#' @description
#'
#' `groupClosest` groups values in `x` for which the difference is smaller than
#' `maxDiff`. As a result, the mean value between the groups will always be
#' larger than `maxDiff`. Values which would be assigned to more than one
#' group, are assigned to the one with the smallest difference to the group
#' mean value.
#'
#' In detail, from the sorted `x`, the function starts from the smallest value
#' defining the first group as the one containing all
#' values in `x` with a difference to this first value which is `<= maxDiff`.
#' The next group is the defined based on the next larger value not being part
#' of the first group and includes all values with a difference `<= maxDiff` to
#' this value. For values fulfilling this criteria but being already part of
#' a previous group, the differences to the mean value of the current group
#' and to the mean of previous groups are compared and values are assigned to
#' the group to which they have the smallest difference.
#'
#' Example: values `1.1, 1.9, 2.2` should be grouped with a `maxDiff = 1`. The
#' first group is defined to include all values for which the difference to the
#' first value (`1.1`) is smaller `maxDiff`. Thus, the first group is defined
#' to contain values `1.1 and 1.9`. Then the next group is defined based on the
#' next larger value not part of any group, `2.2`. This group contains values
#' `1.9` and `2.2` with the value `1.9` being already assigned to the first
#' group. The difference between this value `1.9` and the mean of the
#' current group (`mean(c(1.9, 2.2)`) is then compared to the difference of
#' `1.9` to the mean value of the group `1.9` is already part of
#' (which is `mean(c(1.1, 1.9))`). Since the difference to the second group is
#' smaller, `1.9` is removed from the first group and assigned to the second
#' one.
#'
#' @note
#'
#' This grouping approach, in contrast to [xcms::groupOverlaps()], creates
#' smaller groups and values might not be included into the same group even
#' if their difference is smaller than `maxDiff` (see examples below).
#' 
#' @param x `numeric` of values that should be grouped.
#'
#' @param maxDiff `numeric(1)` defining the threshold for difference between
#'     values in `x` to be grouped into the same group.
#'
#' @return `integer` with the group assignment (values grouped together have
#'     the same value).
#'
#' @author Johannes Rainer
#' 
#' @family grouping operations
#'
#' @export
#' 
#' @examples
#'
#' ## The example described above
#' x <- c(1.1, 1.9, 2.2)
#' groupClosest(x)
#' 
#' x <- c(1.1, 1.5, 1.7, 2.3, 2.7, 4.3, 4.4, 4.9, 5.2, 5.4, 5.8, 6, 7, 9, 9.5, 15)
#'
#' groupClosest(x)
#' ## value 5.2 was initially grouped with 4.3 (because their difference is
#' ## smaller 1, but then re-grouped together with 5.4 because the difference
#' ## between 5.4 (the next value outside the group of 4.3) and 5.2 is smaller
#' ## than its difference to the mean value of the group for value 4.3
#'
#' ## Example for a case in which values are NOT grouped into the same group
#' ## even if the difference between them is <= maxDiff
#' a <- c(4.9, 5.2, 5.4)
#' groupClosest(a, maxDiff = 0.3)
groupClosest <- function(x, maxDiff = 1) {
    if (is.unsorted(x)) {
        idx <- order(x)
        x <- x[idx]
    } else idx <- integer()
    x_len <- length(x)
    x_groups <- rep(NA_integer_, x_len)
    i <- 1
    group_id <- 1
    while(any(is.na(x_groups))) {
        grp <- which(abs(x - x[i]) <= maxDiff)
        ## Check if they are already part of a previous group
        not_in_prev_grp <- is.na(x_groups[grp])
        in_prev_grp <- grp[!not_in_prev_grp]
        ## grp <- grp[not_in_prev_grp]
        if (length(in_prev_grp)) {
            ## compare difference to current x[i] to mean of previous group(s)
            ## they are part of and assign them to the group with the closest
            ## difference
            ## i_diff <- abs(x[in_prev_grp] - x[i])
            ## Compare to the average of the current group.
            i_diff <- abs(x[in_prev_grp] - mean(x[grp]))
            prev_grp <- x_groups[in_prev_grp]
            to_rem <- rep(FALSE, length(in_prev_grp))
            for (j in unique(prev_grp)) {
                j_diff <- abs(x[in_prev_grp] - mean(x[which(x_groups == j)]))
                to_rem <- to_rem | j_diff < i_diff
            }
            grp <- c(in_prev_grp[!to_rem], grp[not_in_prev_grp])
        }
        x_groups[grp] <- group_id
        group_id <- group_id + 1
        i <- which.max(is.na(x_groups))
    }
    x_groups[idx] <- x_groups
    x_groups
}
