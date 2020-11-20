### Function rsq_calc ###
#' Function  rsq_calc
#'
#' Takes a Log2CPM numeric matrix and MArrayLM fit object from lmFit and calculates R-squared for each gene fit.
#'
#' @author John Thompson, \email{jrt@@thompsonclan.org}
#' @references https://support.bioconductor.org/p/17431/  (Mark Robinson)
#' @keywords RNA-Seq; R-squared
#'
#' @param x A normalized log2cpm matrix
#' @param fit A MArrayLM object from lmFit
#'
#' @return A vector of R-squared values for each gene fit.
#'
#' @examples
#'
#'   rsq  = rsq_calc (log2cpm, fitObject)
#'
#' @importFrom assertthat assert_that
#' @importFrom Matrix rowSums
#' @importFrom stringr str_c
#'
#' @export
rsq_calc <- function(x, fit)
{

  assertthat::assert_that(class(x)[[1]] %in% c("data.frame", "matrix"),
                          "MArrayLM" %in% class(fit),
                          is.numeric(as.matrix(x)))

  sst <- Matrix::rowSums(x^2)
  ssr <- sst - fit$df.residual*fit$sigma^2
  return(ssr)
}



