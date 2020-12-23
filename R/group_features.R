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
#' their correlation with each other.
#'
#' Two types of groupings are available:
#' 
#' - `inclusive = FALSE` (the default): the algorithm creates small groups of
#'   highly correlated members, all of which have a correlation with each other
#'   that are `>= threshold`. Note that with this algorithm, rows in `x` could
#'   still have a correlation `>= threshold` with one or more elements of a
#'   group they are not part of. See notes below for more information.
#' - `inclusive = TRUE`: the algorithm creates large groups containing rows that
#'   have a correlation `>= threshold` with at least one element of that group.
#'   For example, if row 1 and 3 have a correlation above the threshold and
#'   rows 3 and 5 too (but correlation between 1 and 5 is below the threshold)
#'   all 3 are grouped into the same group (i.e. rows 1, 3 **and** 5).
#'
#' Note that with parameter `f` it is also possible to pre-define groups of
#' rows that should be further sub-grouped based on correlation with each other.
#' In other words, if `f` is provided, correlations are calculated only between
#' rows with the same value in `f` and hence these pre-defined groups of rows
#' are further sub-grouped based on pairwise correlation. The returned `factor`
#' is then `f` with the additional subgroup appended (and separated with a
#' `"."`). See examples below.
#'
#' @note
#'
#' Implementation note of the grouping algorithm:
#'
#' - all correlations between rows in `x` which are `>= threshold` are
#'   identified and sorted decreasingly.
#' - starting with the pair with the highest correlation groups are defined:
#' - if none of the two is in a group, both are put into the same new group.
#' - if one of the two is already in a group, the other is put into the same
#'   group if **all** correlations of it to that group are `>= threshold`
#'   (and are not `NA`).
#' - if both are already in the same group nothing is done.
#' - if both are in different groups: an element is put into the group of the
#'   other if a) all correlations of it to members of the other's group 
#'   are not `NA` and `>= threshold` **and** b) the average correlation to the
#'   other group is larger than the average correlation to its own group.
#'
#' This ensures that groups are defined in which all elements have a correlation
#' `>= threshold` with each other and the correlation between members of the
#' same group is maximized.
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
#' @param inclusive `logical(1)` whether a version of the grouping algorithm
#'     should be used that leads to larger, more loosely correlated, groups. The
#'     default is `inclusive = FALSE`. See description for more information.
#' 
#' @return `factor` with same length than `nrow(x)` with the group each row
#'     is assigned to.
#'
#' @author Johannes Rainer
#'
#' @importFrom stats cor
#'
#' @importFrom MsFeatures groupSimilarityMatrix
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
                               threshold = 0.9, f = NULL, inclusive = FALSE) {
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
                cors <- cor(t(x[idx, ]), method = method, use = use)
                if (inclusive) {
                    ## Ensure diagonal matrix is TRUE so that even if some
                    ## features have a correlation value of NA they are not
                    ## dropped
                    cors <- cors >= threshold
                    cors[cbind(1:idxl, 1:idxl)] <- TRUE
                    fids <- .index_list_to_factor(.group_logic_matrix(cors))
                } else
                    fids <- groupSimilarityMatrix(cors, threshold = threshold)
                fnew[idx] <- paste0(fg, ".", MsFeatures:::.format_id(fids))
            } else
                fnew[idx] <- paste0(fg, ".001")
        }
        as.factor(fnew)
    } else {
        cors <- cor(t(x), method = method, use = use)
        if (inclusive)
            .index_list_to_factor(.group_logic_matrix(cors >= threshold))
        else
            as.factor(groupSimilarityMatrix(cors, threshold = threshold))
    }
}

#' @title Group EICs based on their correlation
#'
#' @description
#'
#' `groupEicCorrelation` groups (extracted ion) chromatograms (EICs) based on
#' their correlation with each other. If this correlation is `>=` than the
#' provided `threshold` they are grouped.
#'
#' If `x` is a [MChromatograms()] object with more than one column (sample),
#' pairwise correlations between EICs are first calculated for each column
#' (sample) of `x` separately and subsequently aggregated across samples using
#' `aggregationFun`. If `x` is a `MChromatograms` with 4 rows (EICs) and 3
#' columns (samples), pairwise correlations are first calculated between all
#' 4 EICs in each of the 3 columns resulting in 3 correlation matrices (of
#' dimension 4x4). These correlation matrices are combined into a single matrix
#' by combining the 3 correlation values per comparison with
#' `aggregationFun`. By default the mean of the correlation value between e.g.
#' EIC 1 and EIC 2 in each of the 3 columns is used as the final correlation
#' value. Similar to the one-column case EICs are grouped if their (aggregated)
#' correlation coefficient is larger than `threshold`.
#'
#' Two types of groupings are available:
#' 
#' - `inclusive = FALSE` (the default): the algorithm creates small groups of
#'   highly correlated members, all of which have a correlation with each other
#'   that are `>= threshold`. Note that with this algorithm, rows in `x` could
#'   still have a correlation `>= threshold` with one or more elements of a
#'   group they are not part of. See notes below for more information.
#' - `inclusive = TRUE`: the algorithm creates large groups containing rows that
#'   have a correlation `>= threshold` with at least one element of that group.
#'   For example, if row 1 and 3 have a correlation above the threshold and
#'   rows 3 and 5 too (but correlation between 1 and 5 is below the threshold)
#'   all 3 are grouped into the same group (i.e. rows 1, 3 **and** 5).
#'
#' For more information see [groupByCorrelation()].
#'
#' Note that it might be useful to set `tolerance = 0` if chromatograms from
#' the **same** sample are compared. This forces retention times of the compared
#' chromatograms' intensities to be identical.
#' 
#' @param x [MChromatograms()] object of `list` of [Chromatogram()] objects.
#'
#' @param aggregationFun `function` to combine the correlation values between
#'     pairs of EICs across samples (columns). See description for details.
#'
#' @param threshold `numeric(1)` with the threshold for correlation above which
#'     EICs are grouped together.
#'
#' @param align `character(1)` defining the method how chromatograms should be
#'     aligned prior correlation. Defaults to `align = "closest"`. See
#'     [alignRt()] for more details.
#'
#' @param inclusive `logical(1)` defining the grouping approach. With
#'     `inclusive = FALSE` (the default) small groups of highly correlated
#'     features are created using the [groupSimilarityMatrix()] function. With
#'     `inclusive = TRUE` groups are created with features that have at least
#'     one correlation with any other member of the group which is higher than
#'     `threshold`.
#' 
#' @param ... parameters for the [correlate()] function for [MChromatograms()]
#'     objects, such as `tolerance` to allow specifying the maximal acceptable
#'     difference in retention times between objects. See also [alignRt()] for
#'     more information.
#'
#' @return `factor` same length as `nrow(x)` (if `x` is a `MChromatograms`
#'     object) or `length(x)` (if `x` is a `list`) with the group each EIC
#'     is assigned to.
#' 
#' @importMethodsFrom xcms correlate
#' 
#' @importFrom MSnbase MChromatograms
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
#' chrs <- MChromatograms(list(chr1, chr2, chr3))
#'
#' groupEicCorrelation(chrs)
#'
#' ## With a MChromatograms with two columns, use the maximal correlation
#' ## coefficient found in each of the columns
#' chrs <- MChromatograms(list(chr1, chr2, chr3, chr1, chr2, chr3), ncol = 2)
#' groupEicCorrelation(chrs, aggregationFun = max)
groupEicCorrelation <- function(x, aggregationFun = mean,
                                threshold = 0.8, align = "closest",
                                inclusive = FALSE, ...) {
    nr <- nrow(x)
    nc <- ncol(x)
    res <- array(NA_real_, dim = c(nr, nr, nc))
    ## For performance issues it would also be possible to run with full = FALSE
    for (i in seq_len(nc))
        res[, , i] <- correlate(x[, i], align = align, ...)
    res <- apply(res, c(1, 2), aggregationFun, na.rm = TRUE)
    ## Ensure diagonal is always TRUE to not drop any features!
    res[cbind(1:nr, 1:nr)] <- 1
    if (inclusive) {
        res <- res > threshold
        .index_list_to_factor(.group_logic_matrix(res))
    }
    else as.factor(groupSimilarityMatrix(res, threshold = threshold))
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

