#' @include hidden_aliases.R
NULL

#' @rdname hidden_aliases
#'
#' @exportMethod groupFeatures
setGeneric("groupFeatures", function(object, param, ...)
    standardGeneric("groupFeatures"))

#' @title Group features based on approximate retention times
#'
#' @name groupFeatures-approximate-rtime
#'
#' @description
#'
#' Group features based on similar retention time. This method is supposed to be
#' used as an initial *crude* grouping of features based on the median retention
#' time of all their chromatographic peaks. All features with a difference in
#' their retention time which is `<=` parameter `diffRt` of the parameter object
#' are grouped together. Two different grouping methods are available:
#'
#' - `method = "greedy"`: this approach consecutively groups elements together
#'   if their difference in retention time is smaller than `diffRt`. If two
#'   features are grouped into one group, also all other features with a
#'   retention time within the defined window to any of the two features are
#'   also included into the feature group. This grouping is recursively
#'   expanded which can lead, depending on `diffRt` to very large feature
#'   groups spanning a large retention time window.
#'
#' - `method = "groupClosest"`: this approach uses the [groupClosest()] function
#'   that groups values together if their difference is smaller than `diffRt`.
#'   If the difference of a feature to more than one group is smaller `diffRt`
#'   it is assigned to the group to which its retention time is closest (most
#'   similar) to the mean retention time of that group. This leads to smaller
#'   group sizes. See [groupClosest()] for details and examples.
#'
#' @param diffRt `numeric(1)` defining the retention time window within which
#'     features should be grouped. All features with a rtime difference
#'     smaller or equal than `diffRt` are grouped.
#'
#' @param method `character(1)` defining which grouping approach should be
#'     taken. Allowed values are `method = "groupClosest"` (the default) and
#'     `method = "greedy"`. See description for details.
#' 
#' @param msLevel `integer(1)` defining the MS level on which the features
#'     should be grouped.
#' 
#' @param object [XCMSnExp()] object containing also correspondence results.
#' 
#' @param param `SimilarRtimeParam` object with the settings for the method.
#'
#' @return input `XCMSnExp` with feature groups added (i.e. in column
#'     `"feature_group"` of its `featureDefinitions` data frame. 
#'
#' @family feature grouping methods
#' 
#' @rdname groupFeatures-approximate-rtime
#'
#' @importClassesFrom xcms Param
#' 
#' @exportClass SimilarRtimeParam
#'
#' @author Johannes Rainer
NULL

setClass("SimilarRtimeParam",
         slots = c(diffRt = "numeric",
                   method = "character"),
         contains = "Param",
         prototype = prototype(
             diffRt = 1,
             method = "groupClosest"
         ),
         validity = function(object) {
             msg <- NULL
             if (length(object@diffRt) != 1 || object@diffRt < 0)
                 msg <- c("'diffRt' has to be a positive numeric of length 1")
             if (length(object@method) != 1 || !object@method %in%
                 c("greedy", "groupClosest"))
                 msg <- c(
                     msg,
                     "'method' should be either \"groupClosest\" or \"greedy\"")
             msg
         })

#' @rdname groupFeatures-approximate-rtime
#'
#' @export
SimilarRtimeParam <- function(diffRt = 1, method = c("groupClosest", "greedy")) {
    method <- match.arg(method)
    new("SimilarRtimeParam", diffRt = diffRt, method = method)
}

#' @rdname groupFeatures-approximate-rtime
#'
#' @importMethodsFrom xcms hasFeatures featureDefinitions featureDefinitions<-
#'
#' @importFrom MsCoreUtils group
#' 
#' @importClassesFrom xcms XCMSnExp
setMethod(
    "groupFeatures",
    signature(object = "XCMSnExp", param = "SimilarRtimeParam"),
    function(object, param, msLevel = 1L) {
        if (!hasFeatures(object))
            stop("No feature definitions present. Please run ",
                 "first 'groupChromPeaks'")
        if (length(msLevel) > 1)
            stop("Currently only grouping of features from a single MS level",
                 " is supported.")
        is_msLevel <- featureDefinitions(object)$ms_level == msLevel
        if (any(colnames(featureDefinitions(object)) == "feature_group")) {
            f <- featureDefinitions(object)$feature_group
            if (!all(is.na(f[is_msLevel])))
                warning("Found existing feature groups. ",
                        "These will be over-written ")
        } else
            f <- rep(NA_character_, nrow(featureDefinitions(object)))
        if (param@method == "greedy")
            f[is_msLevel] <- paste0(
                "FG.", group(featureDefinitions(object)$rtmed[is_msLevel],
                             tolerance = param@diffRt))
        if (param@method == "groupClosest")
            f[is_msLevel] <- paste0(
                "FG.", groupClosest(featureDefinitions(object)$rtmed[is_msLevel],
                                    maxDiff = param@diffRt))
        featureDefinitions(object)$feature_group <- f
        object
    })

#' @title Group features based on correlation of extracted ion chromatograms
#'
#' @name groupFeatures-eic-correlation
#'
#' @description
#'
#' Group features based on correlation of their extracted ion chromatograms.
#' This correlation is performed separately for each sample with the correlation
#' coefficients being aggregated across samples for the final comparison with
#' parameter `threshold` (the 75% quantile of the per-sample correlation values
#' is used for the comparison with `threshold`).
#'
#' This feature grouping should be called **after** an initial feature
#' grouping by retention time (see [SimilarRtimeParam()]). The feature groups
#' defined in columns `"feature_group"` of `featureDefinitions(object)` (for
#' features matching `msLevel`) will be used and refined by this method.
#' Features with a value of `NA` in `featureDefinitions(object)$feature_group`
#' will be skipped/not considered for feature grouping.
#'
#' While being possible to be performed on the full data set without prior
#' feature grouping , this is not suggested for the following reasons: I) the
#' selection of the top `n` samples with the highest signal for the
#' *feature group* will be biased by very abundant compounds as this is
#' performed on the full data set (i.e. the samples with the highest overall
#' intensities are used for correlation of all features) and II) it is
#' computationally much more expensive because a pairwise correlation between
#' all features has to be performed.
#' 
#' It is also suggested to perform the correlation on a subset of samples
#' per feature with the highest intensities of the peaks (for that feature)
#' although it would also be possible to run the correlation on all samples by
#' setting `n` equal to the total number of samples in the data set. EIC
#' correlation should however be performed ideally on samples in which the
#' original compound is highly abundant to avoid correlation of missing values
#' or noisy peak shapes as much as possible.
#'
#' By default also the signal which is outside identified chromatographic peaks
#' is excluded from the correlation  (parameter `clean`).
#' 
#' @param clean `logical(1)` whether the correlation should be performed only
#'     on the signals within the identified chromatographic peaks
#'     (`clean = TRUE`, default) or all the signal from the extracted ion
#'     chromatogram.
#'
#' @param msLevel `integer(1)` defining the MS level on which the features
#'     should be grouped.
#' 
#' @param n `numeric(1)` defining the total number of samples per feature group
#'     on which this correlation should be performed. This value is rounded up
#'     to the next larger integer value (i.e.  
#' 
#' @param object [XCMSnExp()] object containing also correspondence results.
#' 
#' @param param `EicCorrelationParam` object with the settings for the method.
#'
#' @param threshold `numeric(1)` with the minimal required correlation
#'     coefficient to group featues.
#' 
#' @param value `character(1)` defining whether samples should be grouped based
#'     on the sum of the maximal peak intensity (`value = "maxo"`, the default)
#'     or the integrated peak area (`value = "into"`) for a feature.
#' 
#' @return input `XCMSnExp` with feature groups added (i.e. in column
#'     `"feature_group"` of its `featureDefinitions` data frame. 
#'
#' @family feature grouping methods
#' 
#' @rdname groupFeatures-eic-correlation
#'
#' @exportClass EicCorrelationParam
#'
#' @author Johannes Rainer
NULL

setClass("EicCorrelationParam",
         slots = c(threshold = "numeric",
                   n = "numeric",
                   clean = "logical",
                   value = "character"),
         contains = "Param",
         prototype = prototype(
             threshold = 0.9,
             n = 1,
             clean = TRUE,
             value = "maxo"
         ),
         validity = function(object) {
             msg <- NULL
             if (length(object@threshold) != 1 || object@threshold < 0)
                 msg <- "'threshold' has to be a positive numeric of length 1"
             if (length(object@n) != 1 || object@n < 0)
                 msg <- c(msg, "'n' has to be a positive numeric of length 1")
             if (length(object@clean) != 1)
                 msg <- c(msg, "'clean' has to a logical of length 1")
             if (length(object@value) != 1 && !(object@value %in%
                                                c("maxo", "into")))
                 msg <- c(msg, "'value' has to be either \"maxo\" or \"into\"")
             msg
         })

#' @rdname groupFeatures-eic-correlation
#'
#' @export
EicCorrelationParam <- function(threshold = 0.9, n = 1, clean = TRUE,
                                value = c("maxo", "into")) {
    value = match.arg(value)
    new("EicCorrelationParam", threshold = threshold, n = n, clean = clean,
        value = value)
}

#' @rdname groupFeatures-eic-correlation
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @importMethodsFrom MSnbase fileNames
#'
#' @importMethodsFrom xcms filterFile featureValues removeIntensity
#'
#' @importFrom xcms featureChromatograms
setMethod(
    "groupFeatures",
    signature(object = "XCMSnExp", param = "EicCorrelationParam"),
    function(object, param, msLevel = 1L) {
        if (!hasFeatures(object))
            stop("No feature definitions present. Please run ",
                 "first 'groupChromPeaks'")
        if (length(msLevel) > 1)
            stop("Currently only grouping of features from a single MS level",
                 " is supported.")
        n <- ceiling(param@n)
        nc <- length(fileNames(object))
        if (n > nc)
            stop("'n' should be smaller or equal than the number of ",
                 "samples (", nc, ")")
        is_msLevel <- featureDefinitions(object)$ms_level == msLevel
        if (any(colnames(featureDefinitions(object)) == "feature_group")) {
            f <- featureDefinitions(object)$feature_group
            f_new <- as.character(f)
        } else {
            f <- rep("FG.1", nrow(featureDefinitions(object)))
            f_new <- rep(NA_character_, length(f))
        }
        f[!is_msLevel] <- NA
        if (is.factor(f))
            f <- droplevels(f)
        else
            f <- factor(f, levels = unique(f))
        fvals <- featureValues(object, method = "maxint", value = param@value)
        ffun <- function(z, na.rm = TRUE)
            quantile(z, probs = 0.75, na.rm = na.rm)
        pb <- txtProgressBar(min = 0, max  = nrow(fvals), style = 3)
        setTxtProgressBar(pb, 0)
        counter <- 0
        for (fg in levels(f)) {
            idx <- which(f == fg)
            idxl <- length(idx)
            counter <- counter + idxl
            if (idxl > 1) {
                vals <- apply(fvals[idx, ], MARGIN = 2, sum, na.rm = TRUE)
                sample_idx <- order(vals, decreasing = TRUE)[seq_len(n)]
                obj_sub <- filterFile(object, sample_idx, keepFeatures = TRUE)
                ## Can happen that some of the features are not present in the
                ## subsetted object. Will put them into their own individual
                ## groups later.
                idx_miss <- which(!rownames(fvals)[idx] %in%
                                  rownames(featureDefinitions(obj_sub)))
                if (length(idx_miss)) {
                    idx_miss <- idx[idx_miss]
                    idx <- idx[-idx_miss]
                }
                if (length(idx) > 1) {
                    eics <- featureChromatograms(
                        obj_sub, features = rownames(fvals)[idx], filled = TRUE)
                    if (param@clean)
                        eics <- removeIntensity(eics, which = "outside_chromPeak")
                    res <- groupEicCorrelation(
                        as(eics, "Chromatograms"), aggregationFun = ffun,
                        threshold = param@threshold)
                } else res <- factor(1)
                f_new[idx] <- paste0(fg, ".", res)
                if (length(idx_miss))
                    f_new[idx_miss] <- paste0(
                        fg, ".", seq((length(levels(res)) + 1),
                                     length.out = length(idx_miss)))
            } else
                f_new[idx] <- paste0(fg, ".1")
            setTxtProgressBar(pb, counter)
        }
        setTxtProgressBar(pb, nrow(fvals))
        close(pb)
        featureDefinitions(object)$feature_group <- f_new
        object
    })
