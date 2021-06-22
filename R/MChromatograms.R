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

## split a MChromatograms object.

## take a list of Chromatogram objects of a single column MChromatograms and
## combine the provided Chromatogram objects into a single one.
## combineChromatograms
