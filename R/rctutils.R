# RCT utilities
#

### Function eBayes_autoprop ###
#' eBayes_autoprop
#'
#' Variant of `limma::eBayes()` that sets the `proportion` argument automatically.
#'
#' Limma's [limma::eBayes()] function can calculate log-odds scores
#' (a.k.a. B statistics), which depend on an estimated proportion of
#' differentially expressed genes that the user specifies *a priori*.
#' Limma also provides methods for calculating this proportion via
#' [limma::propTrueNull()], and other packages also provide similar
#' functions. This function combines the two in order to set the
#' proportion automatically.
#'
#' The B statistics produced when running `limma::topTable()` on the
#' value returned by this function should be directly interpretable as
#' estimated log odds, as long as the distribituion of p-values is
#' well-behaved.
#'
#' Note that if the p-value distribution is highly atypical, a
#' proportion that is not between 0 and 1 could be estimated. If this
#' happens, an error will be thrown.
#' @author Ryan Thompson, \email{rct@thompsonclan.org}
#' @param prop.method The method by which to calculate the proportion
#'     of true null hypotheses (i.e. non-differentially expressed
#'     genes). This can be a string to be used as the `method`
#'     argument to [limma::propTrueNull()], or a function that accepts
#'     a vector of p-values and returns a single value between 0 and
#'     1. This argument can only be passed by name.
#' @param ... All other arguments are passed to [limma::eBayes()].
#'
#' @return See [limma::eBayes()].
#'
#' @examples
#'
#' #see ?limma::eBayes
#'
#' @importFrom limma eBayes propTrueNull
#' @importFrom assertthat assert_that
#'
#' @export
eBayes_autoprop <- function(..., prop.method = "lfdr") {
  # req_ns("limma")
  eb <- limma::eBayes(...)
  if (is.function(prop.method)) {
    ptn <- prop.method(eb$p.value)
  } else {
    ptn <- limma::propTrueNull(eb$p.value, method = prop.method)
  }
  assert_that(ptn > 0,ptn < 1)
  limma::eBayes(..., proportion = 1-ptn)
}


### Function get_mds ###
#' get_mds
#'
#' Get a table of MDS values, with proper column names.
#'
#' This runs [limma::plotMDS()], but suppresses the generation of the
#' plot and instead returns the MDS dimensions in a matrix.
#'
#' @author Ryan Thompson, \email{rct@thompsonclan.org}
#' @param x The object to run [limma::plotMDS()] on.
#' @param k The number of MDS dimensions to return. If not specified,
#'     the maximum possible number will be returned.
#' @param ... Additional arguments to [limma::plotMDS()].
#'
#' @return A matrix with `k` columns, and `ncol(x)` rows containing
#'     the MDS dimensions, with each column named "DimN", where N is
#'     the number of that dimension.
#'
#' @examples
#'
#' # TODO Steal from plotMDS
#'
#' @export
get_mds <- function(x, k, ...) {
  req_ns("limma")
  dmat <- limma::plotMDS(x, ..., plot = FALSE)$distance.matrix %>% as.dist
  max_k <- attr(dmat, "Size") - 1
  if (missing(k)) {
    k <- attr(dmat, "Size") - 1
  } else if (k > max_k) {
    warning(glue("Number of requested dimensions ({k}) is greater than the number available ({max_k}). Returning all dimensions."))
    k <- max_k
  }
  mds <- cmdscale(dmat, k = k, eig = TRUE)
  mds$points %>% add_numbered_colnames("Dim")
}


### Function add_numbered_colnames ###
#' add_numbered_colnames
#'
#' Add numbered colnames with a common prefix.
#'
#' @author Ryan Thompson, \email{rct@thompsonclan.org}
#' @param x The variable to add column names to.
#' @param prefix The prefix to use for each column name. Each column's
#'     number will be appended to this prefix to generate the column
#'     name.
#'
#' @return `x`, with its column names set to the column number
#'     appended to `prefix`.
#'
#' @importFrom glue glue
#' @export
add_numbered_colnames <- function(x, prefix = "C") {
  x %>% set_colnames(glue("{prefix}{num}", num = seq(from = 1, length.out = ncol(x))))
}
