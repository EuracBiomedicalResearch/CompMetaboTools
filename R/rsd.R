#' @title Calculate relative standard deviations
#'
#' `rsd` and `rowRsd` are convenience functions to calculate the relative
#' standard deviation (i.e. coefficient of variation) of a numerical vector
#' or for rows of a numerical matrix, respectively.
#'
#' @param x for `rsd` a `numeric`, for `rowRsd` a numerical `matrix`.
#'
#' @param na.rm `logical(1)` whether missing values (`NA`) should be removed
#'     prior to the calculations.
#'
#' @param mad `logical(1)` whether the *Median Absolute Deviation* (MAD) should
#'     be used instead of the standard deviation. This is suggested for
#'     non-gaussian distributed data.
#' 
#' @author Johannes Rainer
#'
#' @md
#'
#' @export
#'
#' @importFrom stats sd
#'
#' @examples
#'
#' a <- c(4.3, 4.5, 3.6, 5.3)
#' rsd(a)
#'
#' A <- rbind(a, a, a)
#' rowRsd(A)
rsd <- function(x, na.rm = TRUE, mad = FALSE) {
    if (mad)
        mad(x, na.rm = na.rm) / abs(median(x, na.rm = na.rm))
    else
        sd(x, na.rm = na.rm) / abs(mean(x, na.rm = na.rm))
}
#' @rdname rsd
#'
#' @export
rowRsd <- function(x, na.rm = TRUE, mad = FALSE)
    apply(x, MARGIN = 1, rsd, na.rm = na.rm, mad = mad)
