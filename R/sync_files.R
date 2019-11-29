#' @title Create local copies of files
#'
#' @description
#'
#' Copy `files` from `remote_path` to `local_path`. Parameter `files` could e.g.
#' be the files defined in a phenodata text file.
#'
#' This functions, along with the related `remove_local_files`, help to
#' manage local copies of e.g. mzML files of an experiment and keep them in sync
#' with a central file storage. A typical use case for these functions is:
#' mzML files are kept in a central storage. To setup or run a specific analysis
#' some of these files (those belonging to the experiment) should be copied to
#' a local copy to enable also off-line analyses. Copying files manually can be
#' cumbersome and error prone, this function thus helps to identify and copy
#' the required files. After an analysis local files might again be removed
#' using the `remove_local_files` function.
#'
#' @param files: `character` with the file names (including relative paths).
#'
#' @param remote_path: `character(1)` with the remote location of the files.
#'
#' @param local_path: `character(1)` local path where to store the files.
#'
#' @author Johannes Rainer
#'
#' @export
#' 
#' @examples
#'
#' ##fls <- c("2018/2018_02/20180203_20000630_NEG.mzML",
#' ##    "2018/2018_02/20180203_20000703_NEG.mzML",
#' ##    "2018/2018_02/20180203_20000735_NEG.mzML")
#'
#' ##sync_files_local(fls, remote_path = "/Users/jo/data",
#' ##    local_path = "/Users/jo/tmp")
sync_files_local <- function(files, remote_path, local_path) {
    if (missing(files) | missing(remote_path) | missing(local_path))
        stop("'files', 'remote_path' and 'local_path' have to be defined")
    files_remote <- paste0(remote_path, "/", files)
    if (!all(file.exists(files_remote)))
        stop("Can not find remote files ",
             paste0(files_remote[!file.exists(files_remote)], collapse = "\n"))
    files_local <- paste0(local_path, "/", files)
    mysub <- !file.exists(files_local)
    files_local <- files_local[mysub]
    files_remote <- files_remote[mysub]
    dirs_local <- unique(dirname(files_local))
    sapply(dirs_local, dir.create, showWarnings = FALSE, recursive = TRUE)
    file.copy(files_remote, files_local, copy.mode = TRUE,
              copy.date = TRUE)
}

#' @rdname sync_files_local
#'
#' @export
remove_local_files <- function(files, local_path) {
    files_remove <- paste0(local_path, "/", files)
    file.remove(files_remove)
}
