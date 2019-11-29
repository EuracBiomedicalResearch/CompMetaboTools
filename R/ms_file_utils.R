#' @title Extract run start time stamp
#'
#' @description
#'
#' Extract the start time stamps from (mzML) files.
#'
#' @param x `character` with the file names (full path) or the mzML files.
#'
#' @param format `character(1)` defining the date/time format of the time
#'     stamp. If `NULL` the time stamp will be returned as a `character`.
#'
#' @param BPPARAM Parallel processing setup. See [bpparam()] for details.
#' 
#' @return `character` (or date/time) with the start time stamps.
#'
#' @author Johannes Rainer
#' 
#' @export
#'
#' @importFrom BiocParallel bpparam
#'
#' @importFrom BiocParallel bplapply
#' 
#' @md
extract_time_stamp <- function(x = character(), format = "%Y-%m-%dT%H:%M:%S",
                               BPPARAM = bpparam()) {
    if (!requireNamespace("mzR", quietly = TRUE))
        stop("The use of this function requires package 'mzR'. Please ",
             "install with 'Biobase::install(\"mzR\")'")
    ts <- bplapply(x, function(z) {
        fl <- mzR::openMSfile(z)
        run_info <- mzR::runInfo(fl)
        mzR::close(fl)
        run_info$startTimeStamp
    }, BPPARAM = BPPARAM)
    if (length(format))
        strptime(unlist(ts), format = format)
    else unlist(ts)
}
