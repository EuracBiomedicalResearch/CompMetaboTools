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
#' are grouped together. If a column `"feature_group"` is found in
#' [xcms::featureDefinitions()] this is further sub-grouped by this method.
#'
#' Two different grouping methods are available:
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
#' @seealso feature-grouping for a general overview.
#'
#' @rdname groupFeatures-approximate-rtime
#'
#' @importClassesFrom xcms Param
#' 
#' @exportClass SimilarRtimeParam
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Performing a quick preprocessing of a test data set.
#' library(faahKO)
#' fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'         system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'         system.file('cdf/WT/wt19.CDF', package = "faahKO"))
#'
#' od <- readMSData(fls, mode = "onDisk")
#' xod <- findChromPeaks(
#'     od, param = CentWaveParam(noise = 10000, snthresh = 40,
#'                               prefilter = c(3, 10000)))
#' pdp <- PeakDensityParam(sampleGroups = c(1, 1, 2))
#' xodg <- groupChromPeaks(xod, param = pdp)
#'
#' ## Group features based on similar retention time (i.e. difference <= 2 seconds)
#' xodg_grp <- groupFeatures(xodg, param = SimilarRtimeParam(diffRt = 2))
#'
#' ## Feature grouping get added to the featureDefinitions in column "feature_group"
#' head(featureDefinitions(xodg_grp)$feature_group)
#'
#' table(featureDefinitions(xodg_grp)$feature_group)
#' length(unique(featureDefinitions(xodg_grp)$feature_group))
#'
#' ## Using the "greedy" method to create larger groups
#' xodg_grp <- groupFeatures(xodg,
#'     param = SimilarRtimeParam(diffRt = 2, method = "greedy"))
#' 
#' length(unique(featureDefinitions(xodg_grp)$feature_group))
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
#' @importClassesFrom xcms XCMSnExp XProcessHistory
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
            f_new <- as.character(f)
        } else {
            f <- rep("FG", nrow(featureDefinitions(object)))
            f_new <- rep(NA_character_, length(f))
        }
        f[!is_msLevel] <- NA
        if (is.factor(f))
            f <- droplevels(f)
        else
            f <- factor(f, levels = unique(f))
        for (fg in levels(f)) {
            idx <- which(f == fg)
            idxl <- length(idx)
            if (idxl > 1) {
                if (param@method == "greedy")
                    fids <- group(featureDefinitions(object)$rtmed[idx],
                                  tolerance = param@diffRt)
                if (param@method == "groupClosest")
                    fids <- groupClosest(
                        featureDefinitions(object)$rtmed[idx],
                        maxDiff = param@diffRt)
                f_new[idx] <- paste0(fg, ".", .format_groups(fids))
            } else
                f_new[idx] <- paste0(fg, ".1")
        }
        featureDefinitions(object)$feature_group <- f_new
        xph <- new("XProcessHistory", param = param, date = date(),
                   type = xcms:::.PROCSTEP.FEATURE.GROUPING,
                   fileIndex = 1:length(fileNames(object)),
                   msLevel = as.integer(msLevel))
        object@.processHistory[[(length(object@.processHistory) + 1)]] <- xph
        validObject(object)
        object
    })

.format_groups <- function(x) {
    digits <- ceiling(log10(length(x) + 1L))
    sprintf(paste0("%0", digits, "d"), as.integer(x))
}

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
#' @param greedy `logical(1)` which grouping algorithm should be used: one that
#'     creates small groups of highly correlated features (`greedy = FALSE`, the
#'     default) or whether features should be grouped that have at least one
#'     correlation with any other member of the group in common
#'     (`greedy = TRUE`). See [groupByCorrelation()] for more information.
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
#' @seealso feature-grouping for a general overview.
#'
#' @rdname groupFeatures-eic-correlation
#'
#' @exportClass EicCorrelationParam
#'
#' @author Johannes Rainer
#'
#' @examples
#' 
#' ## Performing a quick preprocessing of a test data set.
#' library(faahKO)
#' fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'         system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'         system.file('cdf/WT/wt19.CDF', package = "faahKO"))
#'
#' od <- readMSData(fls, mode = "onDisk")
#' xod <- findChromPeaks(
#'     od, param = CentWaveParam(noise = 10000, snthresh = 40,
#'                               prefilter = c(3, 10000)))
#' pdp <- PeakDensityParam(sampleGroups = c(1, 1, 2))
#' xodg <- groupChromPeaks(xod, param = pdp)
#'
#' ## Performing a feature grouping based on EIC correlation on a single
#' ## sample
#' xodg_grp <- groupFeatures(xodg, param = EicCorrelationParam(n = 1))
#'
#' table(featureDefinitions(xodg_grp)$feature_group)
#'
#' ## Usually it is better to perform this correlation on pre-grouped features
#' ## e.g. based on approximate retention time.
#' xodg_grp <- groupFeatures(xodg, param = SimilarRtimeParam(diffRt = 4))
#' xodg_grp <- groupFeatures(xodg_grp, param = EicCorrelationParam(n = 1))
#'
#' table(featureDefinitions(xodg_grp)$feature_group)
NULL

setClass("EicCorrelationParam",
         slots = c(threshold = "numeric",
                   n = "numeric",
                   clean = "logical",
                   value = "character",
                   greedy = "logical"),
         contains = "Param",
         prototype = prototype(
             threshold = 0.9,
             n = 1,
             clean = TRUE,
             value = "maxo",
             greedy = FALSE
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
                                value = c("maxo", "into"), greedy = FALSE) {
    value = match.arg(value)
    new("EicCorrelationParam", threshold = threshold, n = n, clean = clean,
        value = value, greedy = FALSE)
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
            f <- rep("FG", nrow(featureDefinitions(object)))
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
                    tmp <- idx[idx_miss]
                    idx <- idx[-idx_miss]
                    idx_miss <- tmp
                }
                if (length(idx) > 1) {
                    eics <- featureChromatograms(
                        obj_sub, features = rownames(fvals)[idx], filled = TRUE)
                    if (param@clean)
                        eics <- removeIntensity(eics, which = "outside_chromPeak")
                    res <- groupEicCorrelation(
                        as(eics, "Chromatograms"), aggregationFun = ffun,
                        threshold = param@threshold, greedy = param@greedy)
                } else res <- factor(1)
                f_new[idx] <- paste0(fg, ".", .format_groups(res))
                if (length(idx_miss))
                    f_new[idx_miss] <- paste0(
                        fg, ".", .format_groups(
                                     seq((length(levels(res)) + 1),
                                         length.out = length(idx_miss))))
            } else
                f_new[idx] <- paste0(fg, ".1")
            setTxtProgressBar(pb, counter)
        }
        setTxtProgressBar(pb, nrow(fvals))
        close(pb)
        featureDefinitions(object)$feature_group <- f_new
        xph <- new("XProcessHistory", param = param, date = date(),
                   type = xcms:::.PROCSTEP.FEATURE.GROUPING,
                   fileIndex = 1:length(fileNames(object)),
                   msLevel = as.integer(msLevel))
        object@.processHistory[[(length(object@.processHistory) + 1)]] <- xph
        validObject(object)
        object
    })

#' @title Group features based on correlation of abundances across samples
#'
#' @name groupFeatures-abundance-correlation
#'
#' @description
#'
#' This methods groups features based on correlation of abundances (i.e.
#' *feature values*) across samples. Parameter `subset` allows to define a sub
#' set of samples on which this correlation should be correlated. It might for
#' example be better to exclude QC samples from this correlation analysis
#' because values of all features are supposed to be constant in these samples
#' and including these could bias the correlation estimation. Also, it might
#' be better to perform the correlation including also *gap filled* values which
#' is possible with parameter `filled = TRUE` on an `object` on which also
#' [xcms::fillChromPeaks()] was called.
#'
#' @param filled `logical(1)` whether filled-in values should be included in
#'     the correlation analysis. Defaults to `filled = TRUE`.
#' 
#' @param greedy `logical(1)` which grouping algorithm should be used: one that
#'     creates small groups of highly correlated features (`greedy = FALSE`, the
#'     default) or whether features should be grouped that have at least one
#'     correlation with any other member of the group in common
#'     (`greedy = TRUE`). See [groupByCorrelation()] for details.
#' 
#' @param intensity `character(1)` passed to the `featureValues` call. See
#'     [featureValues()] for details. Defaults to `intensity = "maxo"`.
#'
#' @param method `character(1)` passed to the `featureValues` call. See
#'     [featureValues()] for details. Defaults to `method = "maxint"`.
#' 
#' @param msLevel `integer(1)` defining the MS level on which the features
#'     should be grouped.
#' 
#' @param object [XCMSnExp()] object containing also correspondence results.
#' 
#' @param param `AbudanceCorrelationParam` object with the settings for the
#'     method.
#'
#' @param subset `integer` or `logical` defining a subset of samples (at least
#'     2) on which the correlation should be performed.
#'
#' @param transform `function` to be applied to the feature intensities prior
#'     correlation. Defaults to `transform = log2`, i.e. log2 transforms the
#'     data before correlation. To use the values *as-is* use
#'     `transform = identity`.
#' 
#' @param threshold `numeric(1)` with the minimal required correlation
#'     coefficient to group featues.
#' 
#' @param value `character(1)` passed to the `featureValues` call. See
#'     [featureValues()] for details. Defaults to `value = "into"`.
#' 
#' @return input `XCMSnExp` with feature groups added (i.e. in column
#'     `"feature_group"` of its `featureDefinitions` data frame. 
#'
#' @family feature grouping methods
#' 
#' @rdname groupFeatures-abundance-correlation
#'
#' @exportClass AbundanceCorrelationParam
#'
#' @author Johannes Rainer
#'
#' @seealso feature-grouping for a general overview.
#' 
#' @examples
#' 
#' ## Performing a quick preprocessing of a test data set.
#' library(faahKO)
#' fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'         system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'         system.file('cdf/WT/wt19.CDF', package = "faahKO"))
#'
#' od <- readMSData(fls, mode = "onDisk")
#' xod <- findChromPeaks(
#'     od, param = CentWaveParam(noise = 10000, snthresh = 40,
#'                               prefilter = c(3, 10000)))
#' pdp <- PeakDensityParam(sampleGroups = c(1, 1, 2))
#' xodg <- groupChromPeaks(xod, param = pdp)
#'
#' ## Group features based on correlation of feature values (integrated
#' ## peak area) across samples. Note that there are many missing values
#' ## in the feature value which influence grouping of features in the present
#' ## data set.
#' xodg_grp <- groupFeatures(xodg,
#'     param = AbundanceCorrelationParam(threshold = 0.8))
#' table(featureDefinitions(xodg_grp)$feature_group)
#'
#' ## Group based on the maximal peak intensity per feature
#' xodg_grp <- groupFeatures(xodg,
#'     param = AbundanceCorrelationParam(threshold = 0.8, value = "maxo"))
#' table(featureDefinitions(xodg_grp)$feature_group)
NULL

setClass("AbundanceCorrelationParam",
         slots = c(threshold = "numeric",
                   method = "character",
                   value = "character",
                   intensity = "character",
                   filled = "logical",
                   subset = "integer",
                   greedy = "logical",
                   transform = "function"),
         contains = "Param",
         prototype = prototype(
             threshold = 0.9,
             method = "maxint",
             value = "into",
             intensity = "maxo",
             filled = TRUE,
             subset = integer(),
             greedy = FALSE,
             transform = log2
         ),
         validity = function(object) {
             msg <- NULL
             msg
         })

#' @rdname groupFeatures-abundance-correlation
#'
#' @export
AbundanceCorrelationParam <- function(threshold = 0.9, value = "into",
                                      method = c("maxint", "medret", "sum"),
                                      intensity = "maxo", filled = TRUE,
                                      subset = integer(), greedy = FALSE,
                                      transform = log2) {
    method <- match.arg(method)
    if (is.logical(subset))
        subset <- which(subset)
    if (is.numeric(subset))
        subset <- as.integer(subset)
    if (!is.integer(subset))
        stop("'subset' has to be either a logical or an integer vector")
    new("AbundanceCorrelationParam", threshold = threshold, value = value,
        method = method, intensity = intensity, filled = filled,
        subset = subset, greedy = greedy, transform = transform)
}

#' @rdname groupFeatures-abundance-correlation
#'
setMethod(
    "groupFeatures",
    signature(object = "XCMSnExp", param = "AbundanceCorrelationParam"),
    function(object, param, msLevel = 1L) {
        if (!hasFeatures(object))
            stop("No feature definitions present. Please run ",
                 "first 'groupChromPeaks'")
        if (length(msLevel) > 1)
            stop("Currently only grouping of features from a single MS level",
                 " is supported.")
        nc <- length(fileNames(object))
        if (!length(param@subset))
            param@subset <- seq_len(nc)
        if (!all(param@subset %in% seq_len(nc)))
            stop("'subset' has to be between 1 and ", nc)
        if (length(param@subset) < 2)
            stop("Can not calculate correlations for less than 2 samples")
        is_msLevel <- featureDefinitions(object)$ms_level == msLevel
        if (any(colnames(featureDefinitions(object)) == "feature_group")) {
            f <- featureDefinitions(object)$feature_group
            f_new <- as.character(f)
        } else {
            f <- rep("FG", nrow(featureDefinitions(object)))
            f_new <- rep(NA_character_, length(f))
        }
        f[!is_msLevel] <- NA
        if (is.factor(f))
            f <- droplevels(f)
        else
            f <- factor(f, levels = unique(f))
        fvals <- featureValues(
            object, method = param@method, value = param@value,
            intensity = param@intensity,
            filled = param@filled)[, param@subset, drop = FALSE]
        res <- groupByCorrelation(
            param@transform(fvals[is_msLevel, ]), method = "pearson",
            use = "pairwise.complete.obs", threshold = param@threshold,
            f = f[is_msLevel], greedy = param@greedy)
        f_new[is_msLevel] <- as.character(res)
        featureDefinitions(object)$feature_group <- f_new
        xph <- new("XProcessHistory", param = param, date = date(),
                   type = xcms:::.PROCSTEP.FEATURE.GROUPING,
                   fileIndex = 1:length(fileNames(object)),
                   msLevel = as.integer(msLevel))
        object@.processHistory[[(length(object@.processHistory) + 1)]] <- xph
        validObject(object)
        object
    })

#' @title Grouping of features
#'
#' @name feature-grouping
#' 
#' @description
#'
#' After correspondence analysis ([xcms::groupChromPeaks()]) the identified
#' features can be further grouped based on different criteria into
#' *feature groups* which can ideally group features representing signals from
#' ions (adducts, isotopes) of the same compound/metabolite.
#'
#' The available options for the `groupFeatures` method are:
#'
#' - Grouping by similar retention times: [SimilarRtimeParam()].
#'
#' - Grouping by similar feature values across samples:
#'   [AbundanceCorrelationParam()].
#'
#' - Grouping by similar peak shape of extracted ion chromatograms:
#'   [EicCorrelationParam()].
#'
#' An ideal workflow grouping features should sequentially perform the above
#' methods (in the listed order).
#'
#' @param object an [XCMSnExp()] object.
#' 
#' @author Johannes Rainer, Mar Garcia-Aloy, Vinicius Veri Hernandes
#'
#' @seealso [plotFeatureGroups()] for visualization of grouped features.
NULL

#' @rdname feature-grouping
#'
#' @export
featureGroups <- function(object) {
    featureDefinitions(object)$feature_group
}

#' @title Plot feature groups in the m/z-retention time space
#'
#' @description
#'
#' `plotFeatureGroups` visualizes defined feature groups in the m/z by
#' retention time space. Features are indicated by points with features from
#' the same feature group being connected by a line. See [featureGroups()]
#' for details on and options for feature grouping.
#'
#' @param x [XCMSnExp()] object with grouped features (i.e. after calling
#'     [groupFeatures()].
#'
#' @param xlim `numeric(2)` with the lower and upper limit for the x-axis.
#'
#' @param ylim `numeric(2)` with the lower and upper limit for the y-axis.
#'
#' @param xlab `character(1)` with the label for the x-axis.
#'
#' @param ylab `character(1)` with the label for the y-axis.
#'
#' @param pch the plotting character. Defaults to `pch = 4` i.e. plotting
#'     features as crosses. See [par()] for more information.
#'
#' @param col color to be used to draw the features. At present only a single
#'     color is supported.
#'
#' @param type plotting type (see [par()]). Defaults to `type = "o"` which
#'     draws each feature as a point and connecting the features of the same
#'     feature group with a line.
#'
#' @param main `character(1)` with the title of the plot.
#'
#' @param featureGroups optional `character` of feature group IDs to draw only
#'     specified feature group(s). If not provided, all feature groups are
#'     drawn.
#'
#' @export
#'
#' @author Johannes Rainer
plotFeatureGroups <- function(x, xlim = numeric(), ylim = numeric(),
                              xlab = "retention time", ylab = "m/z",
                              pch = 4, col = "#00000060", type = "o",
                              main = "Feature groups",
                              featureGroups = character()) {
    if (!inherits(x, "XCMSnExp"))
        stop("'x' is supposed to be an xcms result object ('XCMSnExp')")
    if (!length(featureGroups(x)))
        stop("No feature groups present. Please run 'groupFeatures' first")
    fts <- factor(featureGroups(x))
    if (!length(featureGroups))
        featureGroups <- levels(fts)
    fts <- fts[fts %in% featureGroups]
    fts <- droplevels(fts)
    if (!length(fts))
        stop("None of the specified feature groups found")
    fdef <- featureDefinitions(x)[featureGroups(x) %in% fts, ]
    rts <- split(fdef$rtmed, fts)
    mzs <- split(fdef$mzmed, fts)
    xy <- cbind(
        x = unlist(lapply(rts, function(z) c(z, NA)), use.names = FALSE),
        y = unlist(lapply(mzs, function(z) c(z, NA)), use.names = FALSE))
    if (length(xlim) != 2)
        xlim <- range(unlist(rts, use.names = FALSE))
    if (length(ylim) != 2)
        ylim <- range(unlist(mzs, use.names = FALSE))
    plot(3, 3, pch = NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
    lines(xy, type = type, col = col, pch = pch)
}
