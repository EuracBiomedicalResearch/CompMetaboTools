#' @title Identify rows with a minimum required proportion of non-missing values
#'
#' @description
#'
#' Helper function to identify rows in `x` which have a proportion of
#' non-missing values (i.e. valid measurements) larger than `prop`. Parameter
#' `f` allows to define groups of columns and hence to perform the test
#' separately in each group. Parameter `condiction` allows to combine these
#' per-group tests to a single `logical` value per row, i.e. the default
#' `condition = any` returns `TRUE` if in at least one group defined by `f` the
#' proportion of non-missing values is `> prop`. With `condition = all` more
#' than `prop` non-missing values have to be present in all groups.
#'
#' @param x `numeric` matrix.
#'
#' @param f optional `factor` defining sample groups (columns) of `x`.
#'
#' @param prop `numeric(1)` defining the minimal required proportion of
#'     non-missing values. Has to be a value between 0 and 1.
#'
#' @param condition optional `function` to be used as condition. With
#'     `condition = any` `TRUE` is reported for all rows in which at least
#'     `prop` values are non-NA within *any* of the sample groups. With
#'     `condition = all` this has to be the case for *all* sample groups.
#'
#' @return `logical` vector, same length than `nrow(x)`, `TRUE` for rows that
#'     pass the criteria.
#'
#' @author Johannes Rainer
#'
#' @export
#' 
#' @examples
#'
#' x <- rbind(
#'     c(NA, 3, 4, 1, 3, NA, 4, NA),
#'     c(4, 2, 3, 4, 5, 5, 2, NA),
#'     c(NA, NA, NA, NA, NA, 3, 4, 5))
#'
#' ## which rows have more than 50% non-missing values
#' moreAreValidThan(x, prop = 0.5)
#'
#' ## Same but with a grouping of columns
#' moreAreValidThan(x, prop = 0.5, f = c(1, 1, 1, 1, 1, 2, 2, 2))
#'
#' ## Same, but require it to be true in all groups
#' moreAreValidThan(x, prop = 0.5, f = c(1, 1, 1, 1, 1, 2, 2, 2), condition = all)
moreAreValidThan <- function(x, f = rep(1, ncol(x)), prop = 0.3,
                             condition = any) {
    if (ncol(x) != length(f))
        stop("length of 'group' has to match number of columns of 'x'")
    apply(x, 1, function(z) {
        condition(vapply(split(z, f),
                         function(y) sum(!is.na(y)) > prop * length(y),
                         logical(1)))
    })
}
