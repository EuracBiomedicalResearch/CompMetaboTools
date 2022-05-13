#' @title Average replicated samples in a SummarizedExperiment
#'
#' @description
#'
#' `averageSE` averages replicated samples (columns) in an
#' `SummarizedExperiment`.
#'
#' @param x `SummarizedExperiment`.
#'
#' @param column `character(1)` defining the column in `colData` that defines
#'     replicates.
#'
#' @param mainAssay optional `character(1)` defining an assay on which per-row
#'     standard deviations across replicates should be calculated. The result
#'     will be stored in an additional assay with the name of `mainAssay` and
#'     suffix `"_sd"`.
#'
#' @return The averages `SummarizedExperiment`.
#'
#' @author Johannes Rainer
#' 
#' @md
#'
#' @export
averageSE <- function(x, column = character(), mainAssay = character()) {
    if (!column %in% colnames(colData(x)))
        stop("Column '", "' not found in 'colData' of 'x'")
    f <- factor(colData(x)[, column], levels = unique(colData(x)[, column]))
    ## new colData: take the first element for each replicate.
    cd <- colData(x)[match(levels(f), f), ]
    rownames(cd) <- cd[, column]
    ## loop over the assays and average them.
    a <- lapply(assays(x), function(z) {
        z <- split.data.frame(t(z), f = f)
        z <- do.call(cbind, lapply(z, colMeans, na.rm = TRUE))
        z[is.na(z)] <- NA
        z
    })
    if (length(mainAssay)) {
        tmp <- split.data.frame(t(assay(x, mainAssay)), f = f)
        tmp <- do.call(cbind, lapply(tmp, function(y) {
            apply(y, MARGIN = 2, FUN = sd, na.rm = TRUE)
        }))
        tmp[is.na(tmp)] <- NA
        a[[paste0(mainAssay, "_sd")]] <- tmp
    }
    SummarizedExperiment(assays = a, rowData = rowData(x),
                         colData = cd, metadata = metadata(x))
}
