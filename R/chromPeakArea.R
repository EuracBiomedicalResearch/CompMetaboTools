#' @title Get MS peak area (m/z and rt range) for chromatographic peaks
#'
#' @description
#'
#' Extract the peak area (defined by the m/z and retention time range) of
#' chromatographic peaks. The chromatographic peaks can be defined by their ID
#' (i.e. rownames of the `chromPeaks` matrix) or by m/z and a retention times.
#'
#' The function identifies first chromatographic peaks in `x` that match the
#' input parameters `mz` and `rt` or `peakId` and then defines the m/z and rt
#' range based on the `"mzmin"`, `"mzmax"`, `"rtmin"` and `"rtmax"` of these.
#' If m/z and rt ranges are provided, only chromatographic peaks are considered
#' that have their apex within the region.
#'
#' As a result a `matrix` is returned with columns `"mzmin"`, `"mzmax"`,
#' `"rtmin"` and `"rtmax"` defining the lower and upper boundaries of the
#' region of the peaks. Each row in this matrix represent the result for one
#' element in the input parameter specifying the peaks of interest (e.g. `mz`,
#' `rt` or `peakId`).
#'
#' @param x `XCMSnExp` object with identified chromatographic peaks.
#'
#' @param mz `numeric` or `matrix` with m/z values to identify chromatographic
#'     peaks from which regions should be extracted. If a `matrix` is provided
#'     it is expeted to have 2 columns with the lower and upper m/z boundary. If
#'     a `numeric` is provided parameter `ppm` is also considered.
#'
#' @param rt `numeric` or `matrix` with retention time values to identify
#'     chromatographic peaks from which the peak areas are supposed to be
#'     defined. If a `matrix` is provided it is expected to be a two column
#'     matrix with the lower and upper rt boundary in which chromatographic
#'     peaks should be found. If a `numeric` is provided, parameter `diffRt`
#'     is considered to define the rt region to search for peaks.
#'
#' @param diffRt `numeric(1)` to define the rt region in which peaks should be
#'     identified if a `numeric` is provided with parameter `rt`. The rt region
#'     is defined as `rt - diffRt` to `rt + diffRt`.
#'
#' @param peakId `character` with the IDs of the chromatographic peaks (i.e.
#'     rownames of `chromPeaks(object)`).
#'
#' @param ppm `numeric(1)` to define the m/z region in which peaks should be
#'     identified if `mz` is a `numeric`. The m/z region is defined as
#'     `mz - ppm(mz, ppm)` to `mz + ppm(mz, ppm)`
#'
#' @return `matrix` with columnd `"mzmin"`, `"mzmax"`, `"rtmin"`, `"rtmax"`
#'     with the m/z and retention time ranges calculated on all chromatographic
#'     peaks in `x` matching `mz` and `rt`. Each row representing the result for
#'     each value in `mz` and `rt`.
#'
#' @importFrom MsCoreUtils ppm
#'
#' @importMethodsFrom xcms hasChromPeaks
#' 
#' @export
#' 
#' @author Johannes Rainer
#'
#' @examples
#' 
#' library(faahKO)
#' fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'          system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'          system.file('cdf/WT/wt19.CDF', package = "faahKO"))
#'
#' od <- readMSData(fls, mode = "onDisk")
#' xod <- findChromPeaks(
#'    od, param = CentWaveParam(noise = 10000, snthresh = 40,
#'                              prefilter = c(3, 10000)))
#'
#' chromPeakArea(xod, rt = c(2500, 2700), diffRt = 100,
#'     mz = cbind(c(400, 300), c(500, 400)))
chromPeakArea <- function(x, mz = numeric(), rt = numeric(), diffRt = 40,
                          ppm = 50, peakId = character()) {
    if (length(peakId))
        stop("Not implemented yet")
    if (!hasChromPeaks(x))
        stop("No chromatographic peaks in 'x'. Please run ",
             "`findChromPeaks' first.")
    if (is.null(dim(mz))) {
        mzp <- ppm(mz, ppm)
        mz <- cbind(mz - mzp, mz + mzp)
    }
    if (is.null(dim(rt)))
        rt <- cbind(rt - diffRt, rt + diffRt)
    l <- nrow(mz)
    if (l != nrow(rt))
        stop("length (nrow) of parameters 'mz' and 'rt' have to match")
    res <- matrix(ncol = 4, nrow = l,
                  dimnames = list(NULL, c("mzmin", "mzmax", "rtmin", "rtmax")))
    for (i in seq_len(l)) {
        pks <- chromPeaks(x, mz = mz[i, , drop = FALSE], ppm = 0,
                          rt = rt[i, , drop = FALSE], type = "apex_within")
        if (nrow(pks))
            res[i, ] <- c(range(pks[, c("mzmin", "mzmax")]),
                          range(pks[, c("rtmin", "rtmax")]))
    }
    res
}
