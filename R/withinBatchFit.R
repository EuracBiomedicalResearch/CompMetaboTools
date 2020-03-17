#' @title Fit a model separately to each batch in an experiment
#'
#' @description
#'
#' `withinBatchFit` fits the model `model` to each row of feature abundances
#' (in the assay matrix specified with `assay`) of `x`. This allows to perform
#' a per-batch-separate estimation of injection order dependent signal drifts.
#'
#' @param x `SummarizedExperiment` with the data on which the model should be
#'     fitted (feature-wise). To estimate the injection order dependent signal
#'     drift on only QC samples, 'x' should contain only the data for the QC
#'     samples of the experiment. It has to contain feature abundances and all
#'     required sample annotations (such as a variable for *batch* and one for
#'     the injection index).
#'
#' @param batch `factor` defining the batches.
#'
#' @param assay `character(1)` specifying the
#'
#' @param model model formula describing the model that should be fitted to
#'     the data within each batch.
#'
#' @param method `character(1)` specifying the method that should be used to
#'     fit the model. Can be either `"lmrob"` for robust regression or `"lm"`
#'     for least squares regression.
#'
#' @param minVals `integer(1)` with the minimum required data values for the
#'     model to be fitted.
#'
#' @param log.transform `logical(1)` whether feature abundances should be
#'     `log2` transformed prior fitting of the model (default is
#'     `log.transform = TRUE`).
#'
#' @return `list` of length equal to the number of batches (levels of `batch`)
#'     each representing the model fits for one batch. `res[[1]]` is thus the
#'     result for the first batch and is a `list` of fitted models, one for
#'     each feature.
#'
#' @author Johannes Rainer
#'
#' @export
withinBatchFit <- function(x, batch = x$batch, assay = "norm",
                           model = y ~ injection_idx, method = "lmrob",
                           minVals = 6, log.transform = TRUE) {
    if (!inherits(x, "SummarizedExperiment"))
        stop("'x' should be a 'SummarizedExperiment'")
    if (length(batch) != ncol(x))
        stop("length of 'batch' has to match the number of columns of 'x'")
    if (!any(assayNames(x) == assay))
        stop("assay '", assay, "' not found in 'x'; 'assay' should be one of ",
             paste(assayNames(x), collapse = ", "))
    if (!is.factor(batch))
        batch <- factor(batch)
    res <- lapply(levels(batch), function(z) {
        x_batch <- x[, batch == z]
        y_vals <- assay(x_batch, assay)
        if (log.transform)
            y_vals <- log2(y_vals)
        xcms:::rowFitModel(model, data = as.data.frame(colData(x_batch)),
                           y = y_vals, method = method, minVals = minVals)
    })
    names(res) <- levels(batch)
    res
}
