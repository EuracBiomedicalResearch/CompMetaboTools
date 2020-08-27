#' @title Extract data from MSnbase objects for use in Spectra
#'
#' @description
#'
#' `extractSpectraData` extracts the spectra data (m/z and intensity values
#' including metadata) from [MSnbase::Spectrum1] or [MSnbase::Spectrum2] objects
#' (or `list` of such objects) and returns these as a `DataFrame` that can
#' be used to create a [Spectra::Spectra] object. This function enables thus
#' to convert data from the *old* `MSnbase` package to the newer `Spectra`
#' package.
#'
#' @param x a `list` of [MSnbase::Spectrum] objects or an object extending
#'     [MSnbase::MSnExp] (such as also the [xcms::XCMSnExp] object) or a
#'     [MSnbase::MSpectra] object.
#'
#' @return [S4Vectors::DataFrame] with the full spectrum data that can be
#'     passed along to the [Spectra::Spectra()] function to create such an
#'     object.
#'
#' @author Johannes Rainer
#'
#' @export
#' 
#' @importFrom S4Vectors DataFrame
#'
#' @importFrom IRanges NumericList
#'
#' @importFrom methods is
#' 
#' @examples
#'
#' ## Read an mzML file with MSnbase
#' fl <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML",
#'     package = "msdata")
#' data <- MSnbase::filterRt(MSnbase::readMSData(fl, mode = "onDisk"),
#'     rt = c(1, 6))
#'
#' ## Extract the data as a DataFrame
#' res <- extractSpectraData(data)
#'
#' ## This can be used as an input for the Spectra constructor of the
#' ## Spectra package.
#' sps <- Spectra::Spectra(res)
#' sps
extractSpectraData <- function(x) {
    if (inherits(x, "MSpectra")) {
        df <- DataFrame(do.call(rbind, lapply(x, MSnbase:::.spectrum_header)))
        df <- cbind(df, MSnbase::mcols(x))
        df$mz <- NumericList(lapply(x, function(z) z@mz))
        df$intensity <- NumericList(lapply(x, function(z) z@intensity))
    } else if (is(x, "list") || inherits(x, "SimpleList")) {
        df <- DataFrame(do.call(rbind, lapply(x, MSnbase:::.spectrum_header)))
        df$mz <- NumericList(lapply(x, function(z) z@mz))
        df$intensity <- NumericList(lapply(x, function(z) z@intensity))
    } else if (inherits(x, "MSnExp")) {
        df <- DataFrame(MSnbase::fData(x))
        df$mz <- NumericList(MSnbase::mz(x))
        df$intensity <- NumericList(MSnbase::intensity(x))
    } else stop("'x' should be either a 'list' of 'Spectrum' objects or an ",
                "object extending 'MSnExp' or 'MSpectra'.")
    df
}
