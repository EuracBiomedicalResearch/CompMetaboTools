#' @title Apply per-batch model adjustment of feature abundances
#'
#' `withinBatchAdjust` adjusts feature abundances within each batch based
#' on the model provided with parameter `models` (which are fitted using the
#' [withinBatchFit()] function.
#'
#' @param x `SummarizedExperiment` with the feature abundances that should be
#'     adjusted.
#'
#' @param models `list` of per-batch feature-wise model fits as returned by
#'     [withinBatchFit()].
#'
#' @param batch `factor` defining the batch assignment of samples in `x`.
#'
#' @param assay `character(1)` with the name of the assay matrix in `x` that
#'     should be adjusted. Defaults to `"norm"`.
#'
#' @param log.transform `logical(1)` defining whether the model fit has been
#'     performed in log scale and the adjustment is also to be performed in
#'     log scale. Note that even if `log.transform` is `TRUE` adjusted feature
#'     abundances are returned in natural scale.
#'
#' @param ... additional arguments to be passed to
#'     `xcms:::applyModelAdjustment`.
#'
#' @return `SummarizedExperiment` (input object `x`) with the feature abundances
#'     of assay `assay` adjusted based on the provided models `models`.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @importMethodsFrom SummarizedExperiment assayNames assay colData assay<-
withinBatchAdjust <- function(x, models, batch = x$batch, assay = "norm",
                              log.transform = TRUE, ...) {
    if (!inherits(x, "SummarizedExperiment"))
        stop("'x' should be a 'SummarizedExperiment'")
    if (length(batch) != ncol(x))
        stop("length of 'batch' has to match the number of columns of 'x'")
    if (!any(assayNames(x) == assay))
        stop("assay '", assay, "' not found in 'x'; 'assay' should be one of ",
             paste(assayNames(x), collapse = ", "))
    if (!is.factor(batch))
        batch <- factor(batch)
    tmp <- assay(x, assay)
    if (log.transform)
        tmp <- log2(tmp)
    for (bt in levels(batch)) {
        b <- batch == bt
        tmp[, b] <- xcms:::applyModelAdjustment(
                               y = tmp[, b, drop = FALSE],
                               lmod = models[[bt]],
                               data = as.data.frame(colData(x))[b, ],
                               ...
                           )
    }
    if (log.transform)
        tmp <- 2^tmp
    assay(x, assay) <- tmp
    x
}
