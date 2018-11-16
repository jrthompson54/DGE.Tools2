### Function lowIntFilter ###
#' Function  lowIntFilter
#'
#' Takes a DGEobj or counts matrix as input and applies a combination of low
#' intensity filters. zFPKM, FPK and raw Count filters are supported.  A gene
#' must pass all filters.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq; counts; low intensity
#'
#' @param x A counts matrix or data.frame OR DGEobj with RNA-Seq data
#'   (required).
#' @param zfpkmThreshold  Genes below this threshold are removed
#'   (Default = -3.0).  Set to -Inf to disable
#'   this filter.
#' @param fpkThreshold Genes below this threshold are removed (Default = 5).
#'   Set to 0 to disable this filter.
#' @param countThreshold Genes below this threshold are removed (Default = 10 ).
#'   Set to 0 to disable this filter.
#' @param sampleFraction The proportion of samples that must meet the thresholds
#'   (Default = 0.5). Range >0 and <=1.
#' @param genelength Vector of genelength for rows of x (required unless x is a
#'   DGEobj)
#'
#' @return Same class as input object with low intensity rows removed
#'
#' @examples
#'
#'   myDgeObj  = lowIntFilter (myDgeObj)
#'
#' @importFrom assertthat assert_that
#' @importFrom DGEobj getItem
#'
#' @export
lowIntFilter <- function(x, zfpkmThreshold=-3.0, fpkThreshold=5, countThreshold=10, sampleFraction=0.5, genelength){

  assertthat::assert_that(class(x)[[1]] %in% c("DGEobj", "data.frame", "matrix"))

  if (class(x)[[1]] != "DGEobj")
    assertthat::assert_that(typeof(x) %in% "integer", "double", "numeric",
                            !missing(genelength))

  #need counts
  if (class(x)[[1]] == "DGEobj")
    counts <- DGEobj::getItem(dgeObj, "counts")
  else
    counts <- x

  #save input class to define output class
  xClass <- class(x)[[1]]

  #Low Intensity Filter
  fracThreshold <- 0.5

  #get genelength if dgeobj
  if (xClass == "DGEobj")
    genelength <-DGEobj::getItem(dgeObj, "geneData")$ExonLength

  #get needed data transformations
  fpkm <- convertCounts(counts, unit="fpkm", geneLength = genelength)
  zfpkm <- zFPKM(as.data.frame(fpkm))
  fpk <- convertCounts(counts, unit="fpk", geneLength = genelength)

  #For data from Xpress, Expect -Inf fpkm values for very short (<150bp genes)  These rows should be removed.
  zFPKM[is.infinite(zFPKM)] <- NA
  idxcc <- complete.cases(zFPKM)  #idx flags rows to keep
  counts <- counts[idxcc,]
  if (xClass == "DGEobj") x <- x[idxcc,]

  #keep zfpkm >= -3 in fracThreshold of samples
  idxzfpkm <- zfpkm >= zfpkmThreshold
  frac <- rowSums(idxzfpkm) / ncol(idxzfpkm)
  idx <- frac >= fracThreshold
  counts <- counts[idx,]
  if (xClass == "DGEobj") x <- x[idx,]

  #keep FPK >= 5 in fracThreshold of samples
  idxfpk <- fpk >= fpkThreshold

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



