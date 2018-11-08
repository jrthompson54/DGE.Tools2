### Function lowIntFilter ###
#' Function  lowIntFilter
#'
#' Take a named list of dataframes where each dataframe has the same
#' column names (e.g. a list of topTable dataframes). Extract
#' the named column from each dataframe and return a matrix.
#'
#' The common use case for this is to provide a list of topTable
#' data frames and extract one column from each file to create
#' a matrix of LogRatios or Pvalues.
#'
#' Technically, this should work as long as the requested colName is present
#' in each dataframe.  The default robust = TRUE should be used unless you
#' are absolutely certain each dataframe in the input list has the same row count
#' and row order.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq; counts; low intensity
#'
#' @param x A counts matrix or data.frame OR DGEobj with RNA-Seq data.
#' @param zThreshold Default = -3.0 Genes below this threshold are considered not detected.
#' @param countThreshold Default = 10 Minimum read counts to keep a gene
#' @param sampleFraction The proportion of samples that must meet the thresholds  Default = 0.5
#' @param genelength Vector of genelength for rows of x (required unless x is a DGEobj)
#'
#' @return Same class as input object with low intensity rows removed
#'
#' @examples
#'
#'   myDgeObj  = lowIntFilter (myDgeObj)
#'
#' @importFrom  assertthat assert_that
#'
#' @export
lowIntFilter <- function(x, xThreshold=-3.0, countThreshold=10, sampleFraction=0.5, genelength){

  assertthat::assert_that(class(x)[[1]] %in% c("DGEobj", "data.frame", "matrix"))

  if (class(x)[[1]] != "DGEobj")
    assertthat::assert_that(typeof(x) %in% "integer", "double", "numeric",
                            !missing(genelength))

  #need counts
  if (class(x)[[1]] == "DGEobj")
    counts <- getItem(dgeObj, "counts")
  else
    counts <- x

  #save input class to define output class
  xClass <- class(x)[[1]]

  #Low Intensity Filter
  fracThreshold <- 0.5

  #low expression filter
  counts <- getItem(dgeObj, "counts")
  if (xClass == "DGEobj")
    genelength <-getItem(dgeObj, "geneData")$ExonLength
  fpkm <- convertCounts(counts, unit="fpkm", geneLength = genelength)
  zfpkm <- zFPKM(as.data.frame(fpkm))
  #For data from Xpress, Expect -Inf values very short (<150bp genes)  These rows should be removed.
  zFPKM[is.infinite(zFPKM)] <- NA
  idxcc <- complete.cases(zFPKM)  #idx flags rows to keep
  counts <- counts[idxcc,]
  if (xClass == "DGEobj") x <- x[idxcc,]

  #keep zfpkm >= -3 in fracThreshold of samples
  idxzfpkm <- zfpkm >= -3.0
  frac <- rowSums(idxzfpkm) / ncol(idxzfpkm)
  idx <- frac >= fracThreshold
  counts <- counts[idx,]
  if (xClass == "DGEobj") x <- x[idx,]

  #overlay a mincount filter
  counts <- getItem(dgeObj, "counts")
  # genelength <-getItem(dgeObj, "geneData")$ExonLength
  idxmin <- counts >= countThreshold
  frac <- rowSums(idxmin)/ncol(idxmin)
  idx <- frac >= fracThreshold
  counts <- counts[idx,]
  if (xClass == "DGEobj") x <- x[idx,]

  if (xClass == "DGEobj") return(x) else return(counts)

}



