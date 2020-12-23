## I can rbind MChromatograms, get then a `matrix` back. which would be OK
## in this case. But note that this only works if we have the same number of
## columns!

#' @noRd
#'
#' @description
#'
#' rbind MChromatograms assuming that the same columns are present!
#'
#' @importMethodsFrom Biobase AnnotatedDataFrame fData
#'
#' @importFrom MsCoreUtils rbindFill
.bind_rows_chromatograms <- function(...) {
    lst <- list(...)
    if (length(lst) == 1L)
        lst <- lst[[1L]]
    if (inherits(lst, "MChromatograms"))
        return(lst)
    ## If that fails we're in trouble
    dta <- do.call(rbind, lapply(lst, function(z) z@.Data))
    pd <- lst[[1]]@phenoData
    fd <- AnnotatedDataFrame(do.call(rbindFill, lapply(lst, fData)))
    if (nrow(fd) == 0)
        fd <- AnnotatedDataFrame(data.frame(matrix(ncol = 0, nrow = nrow(dta))))
    rownames(dta) <- rownames(fd)
    colnames(dta) <- rownames(pd)
    res <- new("MChromatograms")
    res@.Data <- dta
    res@phenoData <- pd
    res@featureData <- fd
    validObject(res)
    res
}

#' @title Combine MChromatograms objects
#'
#' @description
#'
#' [MChromatograms()] objects can be concatenated (row-wise) with the `c`
#' function resulting in a `MChromatograms` object with the same number of
#' samples (columns).
#'
#' @param x `MChromatograms` object.
#'
#' @param ... `MChromatograms` object that should be appended to `x`.
#' 
#' @note
#'
#' The `c` function for `XChromatograms` objects (defined in the `xcms` package)
#' is not supported.
#' 
#' @return `MChromatograms` object.
#'
#' @rdname c-MChromatograms
#' 
#' @author Johannes Rainer
#'
#' @export
setMethod("c", "MChromatograms", function(x, ...) {
    .bind_rows_chromatograms(unname(c(list(x), list(...))))
})

#' @rdname c-MChromatograms
#'
#' @export
setMethod("c", "XChromatograms", function(x, ...) {
    stop("Concatenation of 'XChromatogram' objects is currently not supported.")
})

#' @title Overlay plots from multiple EICs
#'
#' @description
#'
#' `plotOverlay` draws chromatographic peak data from multiple (different)
#' extracted ion chromatograms (EICs) into the same plot. This allows to
#' directly compare the peak shape of these EICs in the same sample. In
#' contrast to the `plot` function for [MChromatograms()] object, which draws
#' the data from the same EIC across multiple samples in the same plot, this
#' function draws the different EICs from the same sample into the same plot.
#'
#' @param x [MChromatograms()] object with the EICs to draw.
#'
#' @param col color representation. Should be either of length 1 (to use the
#'     same color for each EIC) or equal to `nrow(x)` to use a different color
#'     for each EIC (row in `x`).
#'
#' @param type `character(1)` specifying the plot type (see help from base
#'     `plot` function). `type = "l"` (the default) draws lines, `type = "p"`
#'     would draw each intensity as a point etc.
#'
#' @param main optional `character` with the title for the plot of each sample.
#'     Can be either of length 1 or equal `ncol(x)`.
#'
#' @param xlab `character(1)` with the label for the x-axis.
#'
#' @param ylab `character(1)` with the label for the y-axis.
#' 
#' @param ... additional parameters to be passed down to the `points` function.
#'
#' @author Johannes Rainer
#'
#' @importFrom graphics par
#'
#' @export
plotOverlay <- function(x, col = "#00000060", type = "l", main = NULL,
                        xlab = "rtime", ylab = "intensity", ...) {
    if (!inherits(x, "MChromatograms"))
        stop("'x' is supposed to be a 'MChromatograms' object")
    nc <- ncol(x)
    nr <- nrow(x)
    if (nc > 1)
        par(mfrow = c(round(sqrt(nc)), ceiling(sqrt(nc))))
    if (length(col) != nr)
        col <- rep(col[1], nr)
    if (is.null(main))
        main <- colnames(x)
    if (length(main) != nc)
        main <- rep(main[1], nc)
    for (i in seq_len(ncol(x))) {
        if (nr == 1)
            data <- list(data.frame(rtime = x[1, i]@rtime,
                                    intensity = x[1, i]@intensity))
        else
            data <- lapply(x[, i], function(z)
                data.frame(rtime = z@rtime, intensity = z@intensity))
        xl <- vapply(data, function(z) range(z$rtime, na.rm = TRUE), numeric(2))
        yl <- vapply(data, function(z) range(z$intensity, na.rm = TRUE),
                     numeric(2))
        plot(3, 3, pch = NA, xlim = range(xl, na.rm = TRUE, finite = TRUE),
             ylim = range(yl, na.rm = TRUE, finite = TRUE), main = main[i],
             xlab = xlab, ylab = ylab, ...)
        for (j in seq_along(data))
            points(data[[j]]$rtime, data[[j]]$intensity, col = col[j],
                   type = type, ...)
    }
}

## split a MChromatograms object.

## take a list of Chromatogram objects of a single column MChromatograms and
## combine the provided Chromatogram objects into a single one.
## combineChromatograms
