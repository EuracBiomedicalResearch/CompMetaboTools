#' @include hidden_aliases.R
NULL

setGeneric("joyPlot", def = function(object, ...)
    standardGeneric("joyPlot"))

#' @title Create a stacked plot of multiple chromatograms
#'
#' @importClassesFrom MSnbase MChromatograms
#'
#' @importClassesFrom xcms XCMSnExp
#' 
#' @description
#'
#' The `joyPlot` function creates a stacked plot of chromatograms of a
#' single sample (parameter `sample`). The name `joyPlot` was chosen as the
#' created plot can, depending on the data, look similar to the cover of the
#' [Unknown Pleasures](https://en.wikipedia.org/wiki/Unknown_Pleasures) album
#' of the band Joy Division.
#'
#' The method is implemented for [XCMSnExp-class] and [MChromatograms-class]
#' objects:
#'
#' - For `object` being an [MChromatograms-class]: the function plots all
#'   chromatograms stacked (overlappingly) on top of each others assuming that
#'   the chromatograms are ordered by m/z (i.e. first row having the smallest,
#'   last row the largest m/z).
#'
#' - For `object` being an [XCMSnExp-class]: the function will plot only
#'   chromatograms for chromatographic peaks present in the specified m/z range
#'   (parameter `mz`) and retention time range (parameter `rt`) in the sample
#'   `sample`.
#'
#' @param object an [MChromatograms-class] or [XCMSnExp-class] object.
#'
#' @param sample `integer(1)` defining the sample from which the plot should be
#'     created.
#'
#' @param mz for `object` being an `XCMSnExp`: `numeric(2)` defining the lower
#'     and upper m/z of the m/z range.
#'
#' @param rt for `object` being an `XCMSnExp`: `numeric(2)` defining the
#'     retention time range.
#'
#' @param sample `integer(1)` defining the sample (column in `x`) from which the
#'     plot should be created.
#' 
#' @param yoffset `numeric(1)` defining the proportion of the y-axis to be used
#'     for the stacked chromatograms. The default value of `yoffset = 0.3` uses
#'     30% of the y axis to place the individual chromatograms (eiher equally
#'     spaced or relative to their m/z). This value is relative to the largest
#'     intensity value of all chromatograms.
#'
#' @param xlim optional `numeric(2)` defining the x-axis limits.
#'
#' @param ylim optional `numeric(2)` defining the y-axis limits.
#'
#' @param col color of the lines used to draw the chromatograms. Can be of
#'     length 1 or equal `nrow(x)`.
#'
#' @param peakBg background color for the chromatograms. Can be of length 1
#'     or equal `nrow(x)`.
#'
#' @param type `character(1)` defining the plot type (see parameter `type` in
#'     base [plot()].
#'
#' @param ylab `character(1)` with the y-axis label.
#'
#' @param xlab `character(1)` with the x-axis label.
#'
#' @param spacing `character(1)` defining whether the chromatograms should be
#'     placed on the y-axis relative to their m/z value (spacing = "relative")
#'     or equally spaced (spacing = "equal").
#'
#' @param legend `logical(1)` whether m/z values should be displayed next to
#'     the chromatograms.
#'
#' @param ... additional arguments to be passed to the `plot` function.
#' 
#' @return besides creating the plot the function returns (invisibly) the
#'     y positions of the individual chromatograms.
#'
#' @author Johannes Rainer
#'
#' @importMethodsFrom xcms filterFile rtime
#' 
#' @name joyPlot
#' 
#' @examples
#'
#' library(xcms)
#'
#' ## Load a test data set with detected peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' one <- filterFile(faahko_sub, 1L)
#'
#' ## Plot stacked chromatograms for m/z slices containing chromatographic peaks
#' joyPlot(one, rt = c(4000, 4200))
#'
#' ## Plot the chromatograms equally spaced and showing also their m/z
#' res <- joyPlot(one, rt = c(4000, 4200), spacing = "equal", legend = TRUE)
#'
#' ## joyPlot returns the y position of the slices/chromatograms
#' res
#'
#' ## Disable transparency of individual peaks
#' joyPlot(one, rt = c(4000, 4200), spacing = "equal", peakBg = "#ffffff")
#'
#' ## To create a joy plot for the full data, not just the slices containing
#' ## chromatographic peaks we have to first extract the respective
#' ## chromatograms from the data. Below we define thus the m/z slices - with
#' ## a rather large m/z range for each (i.e. of 1) in the m/z range 400 - 500
#' mzs <- seq(400, 500, by = 1)
#' chr <- chromatogram(one, mz = cbind(mzs[-length(mzs)], mzs[-1]), include = "none", rt = c(4000, 4200))
#'
#' joyPlot(chr)
#'
#' ## And specyfying a different color for each chromatogram
#' joyPlot(chr, col = heat.colors(nrow(chr)))
#'
#' ## Finally, creating a plot that looks more like the album cover
#' par(bg = "#000000")
#' joyPlot(chr, col = "#ffffff", peakBg = "#000000aa", yoffset = 0.7)
#'
NULL

#' @rdname joyPlot
#'
#' @export
setMethod("joyPlot", "MChromatograms",
          function(object, sample = 1, yoffset = 0.3, xlim = NULL, ylim = NULL,
                   col = "#000000", peakBg = "#ffffffaa", type = "l", ylab = "",
                   xlab = "retention time", spacing = c("relative", "equal"),
                   legend = FALSE, ...) {
              .joy_plot_chromatograms(object, sample = sample,
                                      yoffset = yoffset, xlim = xlim,
                                      ylim = ylim, col = col, peakBg = peakBg,
                                      type = type, ylab = ylab, xlab = xlab,
                                      spacing = spacing, legend = legend, ...)
})

#' @rdname joyPlot
#'
#' @importFrom xcms groupOverlaps
#'
#' @importMethodsFrom xcms chromPeaks intensity mz chromatogram
#' 
#' @export
setMethod("joyPlot", "XCMSnExp",
          function(object, sample = 1L, mz = c(-Inf, Inf), rt = c(-Inf, Inf),
                   yoffset = 0.3, col = "#000000", peakBg = "#ffffffaa",
                   type = "l", ylab = "", xlab = "retention time",
                   spacing = c("relative", "equal"), legend = FALSE, ...) {
              object <- filterFile(object, sample)
              pks <- chromPeaks(object, mz = mz, rt = rt)
              rngs <- groupOverlaps(pks[, "mzmin"], pks[, "mzmax"])
              mzr <- do.call(rbind, lapply(rngs, function(z)
                  range(pks[z, c("mzmin", "mzmax")])))
              .joy_plot_chromatograms(
                  chromatogram(object, mz = mzr, rt = rt, include = "none"),
                  yoffset = yoffset, col = col, peakBg = peakBg, type = type,
                  ylab = ylab, xlab = xlab, spacing = spacing, legend = legend,
                  ...)
          })

#' @description
#' 
#' Note that it is assumed the chromatograms in `x` to be ordered by m/z as they
#' are plotted in reverse order (last chromatogram first).
#'
#' @importFrom graphics axis mtext plot.new plot.window polygon title
#'
#' @importFrom grDevices dev.flush dev.hold
#' @noRd
.joy_plot_chromatograms <- function(x, sample = 1L, yoffset = 0.3, xlim = NULL,
                                    ylim = NULL, col = "#000000",
                                    peakBg = "#ffffffaa",
                                    type = "l", ylab = "",
                                    xlab = "retention time",
                                    spacing = c("relative", "equal"),
                                    legend = FALSE, ...) {
    spacing <- match.arg(spacing)
    x <- x[, sample]
    nrx <- nrow(x)
    if (length(col) == 1) {
        col <- rep(col, nrx)
    } else if (length(col) != nrx) {
        warning("Length of 'col' does not match the number of chromatograms.",
                " Using 'col[1]' for all.", call. = TRUE)
        col <- rep(col[1L], nrx)
    }
    if (length(peakBg) == 1) {
        peakBg <- rep(peakBg, nrx)
    } else if (length(peakBg) != nrx) {
        warning("Length of 'peakBg' does not match the number of ",
                "chromatograms. Using 'peakBg[1]' for all.", call. = TRUE)
        peakBg <- rep(peakBg[1L], nrx)
    }
    mzs <- mz(x)
    if (any(is.na(mzs))) {
        mzs <- seq_len(nrx)
    } else mzs <- rowMeans(mzs)
    ymaxs <- vapply(x, function(z)
        suppressWarnings(max(intensity(z), na.rm = TRUE)), numeric(1))
    if (length(ylim)) {
        ylim <- range(ylim)
        maxint <- ylim[2]
        minint <- ylim[1]
    } else {
        maxint <- max(ymaxs, na.rm = TRUE)
        minint <- 0
    }
    yoffset <- maxint * yoffset
    if (spacing == "equal")
        ypos <- seq(0, yoffset, length.out = nrx)
    else {
        mzs <- mzs - min(mzs)
        ypos <- mzs * yoffset / max(mzs)
    }
    if (length(xlim))
        xlim <- range(xlim)
    else
        xlim <- range(vapply(x, function(z) range(rtime(z), na.rm = TRUE),
                             numeric(2)))
    dev.hold()
    on.exit(dev.flush())
    plot.new()
    plot.window(xlim = xlim, ylim = c(minint, max(ypos + ymaxs, na.rm = TRUE)))
    axis(side = 1, ...)
    title(xlab = xlab, ylab = ylab, ...)
    ## plot(3, 3, pch = NA, xlim = xlim, ylim = c(0, maxint + max(ypos)),
    ##      yaxt = "n", bty = "n", xlab = xlab, ylab = ylab, ...)
    for (i in rev(seq_len(nrx))) {
        tmp <- x[i, ]
        rts <- rtime(tmp)
        ints <- intensity(tmp)
        nnas <- !is.na(ints)
        rts <- rts[nnas]
        ints <- ints[nnas]
        if (length(ints)) {
            polygon(c(rts[1], rts, rts[length(rts)]),
                    c(ypos[i], ints + ypos[i], ypos[i]) - minint, border = NA,
                    col = peakBg[i])
            points(rtime(tmp), intensity(tmp) + ypos[i] - minint, type = type,
                   col = col[i], ...)
        }
    }
    if (legend) {
        mtext(side = 2, at = ypos, text = format(mzs, digits = 6), las = 2)
    }
    invisible(ypos)
}
