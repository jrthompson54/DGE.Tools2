### Function lowIntFilter ###
#' Function  lowIntFilter
#'
#' Takes a DGEobj or counts matrix as input and applies a combination of low
#' intensity filters. raw count, zFPKM and/or FPK filters are supported.  A gene
#' must pass all active filters.  Not setting a threshold argument inactivates that threshold.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq; counts; low intensity
#'
#' @param x A counts matrix or data.frame OR a DGEobj with RNA-Seq data
#'   (required).
#' @param countThreshold Genes below this threshold are removed (recommend 10).
#'   Set to 0 to disable this filter.
#' @param zfpkmThreshold  Genes below this threshold are removed (-3.0 is recommended)
#' @param fpkThreshold Genes below this threshold are removed (recommend 5).
#' @param sampleFraction The proportion of samples that must meet the thresholds
#'   (Default = 0.5). Range >0 and <=1.
#' @param genelength Vector of genelength for rows of x (required unless x is a
#'   DGEobj)
#'
#' @return Same class as input object with low intensity rows removed
#'
#' @examples
#'   #simple count threshold
#'   myDgeObj  = lowIntFilter (myDgeObj, countThreshold=10)
#'
#'   #count and zFPKM thresholds
#'   myDgeObj  = lowIntFilter (myDgeObj, countThreshold=10, zfpkmThreshold = -3.0)
#'
#' @importFrom assertthat assert_that
#' @importFrom DGEobj getItem
#' @importFrom zFPKM zFPKM
#'
#' @export
lowIntFilter <- function(x, zfpkmThreshold, fpkThreshold, countThreshold, sampleFraction=0.5, genelength){

  assertthat::assert_that(class(x)[[1]] %in% c("DGEobj", "data.frame", "matrix"))

  if (class(x)[[1]] != "DGEobj")
    assertthat::assert_that(typeof(x) %in% "integer", "double", "numeric",
                            !missing(genelength))

  #need counts
  if (class(x)[[1]] == "DGEobj")
    counts <-x$counts
  else
    counts <- x

  #save input class to define output class
  xClass <- class(x)[[1]]

  #Define sampleFraction
  if (missing(sampleFraction)) sampleFraction <- 0.5

  #get genelength if dgeobj
  if (xClass == "DGEobj"  &  missing(genelength))  #user supplied genelength supercedes
    genelength <-x$geneData$ExonLength

  #apply count threshold
  if (!missing(countThreshold)){
    #overlay a mincount filter
    idxmin <- counts >= countThreshold
    frac <- rowSums(idxmin)/ncol(idxmin)
    idx <- frac >= sampleFraction
    counts <- counts[idx,]
    if (xClass == "DGEobj") x <- x[idx,]
  }

  #apply zFPKM threshold
  if (!missing(zfpkmThreshold)){
    fpkm <- convertCounts(counts, unit="fpkm", geneLength = genelength)
    zfpkm <- zFPKM::zFPKM(as.data.frame(fpkm))

    #For data from Xpress, Expect -Inf fpkm values for very short (<150bp genes)  These rows should be removed.
    zFPKM[is.infinite(zFPKM)] <- NA
    idxcc <- complete.cases(zFPKM)  #idxcc flags rows to keep; drops rows with NAs

    #creat index for zfpkm >= zfpkmThreshold in fracThreshold of samples
    idxzfpkm <- zfpkm >= zfpkmThreshold
    frac <- rowSums(idxzfpkm) / ncol(idxzfpkm)
    idx <- frac >= fracThreshold
    counts <- counts[(idx & idxcc),]
    if (xClass == "DGEobj") x <- x[(idx & idxcc),]
  }

  #apply FPK threshold
  if (!missing(fpkThreshold)){
    fpk <- convertCounts(counts, unit="fpk", geneLength = genelength)
    #keep FPK >= fpkThreshold in fracThreshold of samples
    idxfpk <- fpk >= fpkThreshold
    counts <- counts[idxfpk,]
    if (xClass == "DGEobj") x <- x[idxfpk,]
  }

  if (xClass == "DGEobj") return(x) else return(counts)

}



