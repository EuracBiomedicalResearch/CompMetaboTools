#' @title Match features based on their retention time and m/z values
#'
#' @description
#'
#' `matchRtMz` matches elements (features) in `x` with elements in `table`
#' based on similarity of their m/z **and** retention time. With parameter
#' `duplicates = "closest"`, the function return the index of the best match,
#' considering the m/z difference of features that are within the acceptable
#' retention time difference (defined by `rt_tolerance`). With
#' `duplicates = "keep"` the index of all matching rows in `table` are returned.
#'
#' @note
#'
#' The function first finds features in `table` with a difference of retention
#' time which is smaller than `rt_tolerance` and matches these using the
#' [closest()] function.
#'
#' @param x `data.frame`, `matrix` or `DataFrame` with feature definitions (i.e.
#'     m/z and retentin times) that should be matched against features in
#'     `table`.
#'
#' @param table `data.frame`, `matrix` or `DataFrame` with feature definitions
#'     to match features in `x` against.
#'
#' @param nomatch value that should be returned if no match for a feature is
#'     found.
#'
#' @param rt_tolerance `numeric(1)` with the largest acceptable difference in
#'     retention time.
#'
#' @param tolerance `numeric(1)` with a constant acceptable difference of m/z
#'     values for features to be considered matching.
#'
#' @param ppm `numeric(1)` with a m/z-dependent relative acceptable difference
#'     (in parts per million) of m/z values.
#
#' @param duplicates `character(1)` whether the best (`duplicates = "closest",
#'     default) or all matches (`duplicates = "keep"`) shpuld be returned.
#'
#' @param mzcol `character(1)` with the name of the column containing the m/z
#'     ratios.
#'
#' @param rtcol `character(1)` with the name of the column containing the
#'     retention times.
#'
#' @return for `duplicates = "closest"`: `integer` of length equal to `nrow(x)`
#'     with the index of the row in `table` that matches each row in `x`
#'     (e.g. `c(3, 4)` means the first feature in `x` matches with the 3rd
#'     feature in `table`.
#'     For `duplicates = "keep"`: `list` of length equal to `nrow(x)` with
#'     indices of **all** rows in `table` that match each row in `x`.
#'
#' @importFrom MsCoreUtils closest
#'
#' @export
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' x <- data.frame(mz = c(23.4, 45.6, 56.9, 76.5, 76.5, 76.5, 80.1),
#'     rt = c(12, 34, 59, 34, 67, 65, 67))
#'
#' set.seed(123)
#' y <- rbind(x, x)
#' y$mz <- y$mz + rnorm(nrow(y), sd = 0.0002)
#' y$rt[1:nrow(x)] <- x$rt + 2
#' y <- y[order(y$mz), ]
#'
#' matchRtMz(x, y)
#'
#' ## Keeping all matches
#' matchRtMz(x, y, duplicates = "keep")
#'
#' ## Lower ppm
#' matchRtMz(x, y, duplicates = "keep", ppm = 5)
#'
#' ## even lower
#' matchRtMz(x, y, duplicates = "keep", ppm = 2)
#'
#' matchRtMz(x, y, ppm = 0)
matchRtMz <- function(x, table, nomatch = NA_integer_, rt_tolerance = 2,
                          tolerance = 0, ppm = 20,
                          duplicates = c("closest", "keep"),
                          mzcol = "mz", rtcol = "rt") {
    duplicates <- match.arg(duplicates)
    if (length(dim(x)) != 2 || length(dim(table)) != 2)
        stop("'x' and 'table' have to be two data frames")
    if (!all(c(mzcol, rtcol) %in% colnames(x)) ||
        !all(c(mzcol, rtcol) %in% colnames(table)))
        stop("Columns '", mzcol, "' and '", rtcol, "' not found in 'x' and ",
             "'table'")
    mz1 <- x[, mzcol]
    rt1 <- x[, rtcol]
    mz2 <- table[, mzcol]
    rt2 <- table[, rtcol]
    idxl <- vector("list", length = nrow(x))
    for (i in seq_along(idxl)) {
        matches <- which(abs(rt2 - rt1[i]) <= rt_tolerance)
        if (length(matches)) {
            matches <- matches[which(closest(mz2[matches], mz1[i],
                                             tolerance = tolerance,
                                             ppm = ppm,
                                             duplicates = duplicates,
                                             nomatch = nomatch) == 1)]
            if (length(matches))
                idxl[[i]] <- matches
            else idxl[[i]] <- nomatch
        } else idxl[[i]] <- nomatch
    }
    if (duplicates == "closest")
        as.integer(idxl)
    else idxl
}
