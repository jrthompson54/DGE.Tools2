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
#' @param x A counts matrix or data.frame OR a DGEobj with RNA-Seq (counts) data
#'   (required).
#' @param countThreshold Genes below this threshold are removed (recommend 10).
#'   Set to 0 to disable this filter.
#' @param zfpkmThreshold  Genes below this threshold are removed (-3.0 is recommended)
#' @param fpkThreshold Genes below this threshold are removed (recommend 5).
#' @param sampleFraction The proportion of samples that must meet the thresholds
#'   (Default = 0.5). Range >0 and <=1.
#' @param genelength Vector of genelengths for rows of x. Required for FPK and
#'   zFPKM filters, unless x is a DGEobj.  In the DGEobj case, genelength is
#'   retrieved from the DGEobj, unless supplied by the genelength argument.
#' @param verbose Prints some messages about the filtering process.
#'
#' @return Same class as input object with low intensity rows removed
#'
#' @examples
#'
#'   #simple count threshold
#'   myDgeObj  = lowIntFilter (myDgeObj, countThreshold=10)
#'
#'   #count and zFPKM thresholds
#'   myDgeObj  = lowIntFilter (myDgeObj, countThreshold=10, zfpkmThreshold = -3.0)
#'
#' @importFrom assertthat assert_that
#' @importFrom zFPKM zFPKM
#' @importFrom stringr str_c
#'
#' @export
lowIntFilter <- function(x, zfpkmThreshold, fpkThreshold, countThreshold,
                         sampleFraction=0.5, genelength=NULL,
                         verbose=FALSE)
{

  assertthat::assert_that(class(x)[[1]] %in% c("DGEobj", "data.frame", "matrix"))

  if (class(x)[[1]] != "DGEobj")
    assertthat::assert_that(typeof(x) %in% c("integer", "double", "numeric"))

  #need counts
  if (class(x)[[1]] == "DGEobj")
    counts <-x$counts
  else
    counts <- x

  starting_rowcount <- nrow(counts)

  #save input class to define output class
  xClass <- class(x)[[1]]

  #get genelength if dgeobj
  if (xClass == "DGEobj"  &  is.null(genelength))  #user supplied genelength supercedes genelength from DGEobj
    genelength <-x$geneData$ExonLength

  #apply zFPKM threshold
  if (!missing(zfpkmThreshold)){
    assertthat::assert_that(!is.null(genelength),
                            length(genelength) == nrow(counts))
    fpkm <- convertCounts(counts, unit="fpkm", geneLength = genelength) #class(fpkm) = "matrix"

    #Xpress counts produce rows filled with NaNs.  Need to filter these out first before calculating zFPKM.
    idx <- complete.cases(fpkm)
    fpkm <- fpkm[idx,]
    genelength <- genelength[idx]
    nan_genes <- sum(!idx)
    if (nan_genes > 0 & verbose==TRUE)
      message(stringr::str_c(nan_genes, " genes with NaN FPKM values removed."))

    #calculate zfpkm
    zfpkm <-as.matrix(zFPKM::zFPKM(as.data.frame(fpkm)))

    #create index for zfpkm >= zfpkmThreshold in fracThreshold of samples
    idxzfpkm <- zfpkm >= zfpkmThreshold
    frac <- rowSums(idxzfpkm) / ncol(idxzfpkm)
    idx <- frac >= sampleFraction
    counts <- counts[idx,]
    if (xClass == "DGEobj") x <- x[idx,]
    genelength <- genelength[idx]

    if (verbose == TRUE)
      message(stringr::str_c(sum(idx), " of ", starting_rowcount, " genes retained by the zFPKM filter."))
  }

  #apply FPK threshold
  if (!missing(fpkThreshold)){
    assertthat::assert_that(!is.null(genelength),
                            length(genelength) == nrow(counts))
    fpk <- convertCounts(counts, unit="fpk", geneLength = genelength)
    #keep FPK >= fpkThreshold in fracThreshold of samples
    idxfpk <- fpk >= fpkThreshold
    frac <- rowSums(idxfpk)/ncol(idxfpk)
    idx <- frac >= sampleFraction
    counts <- counts[idx,]
    if (xClass == "DGEobj") x <- x[idx,]
    genelength <- genelength[idx]

    if (verbose == TRUE)
      tsmsg(stringr::str_c(sum(idx), " of ", length(idx), " genes retained by the FPK filter."))
  }

  #apply count threshold
  if (!missing(countThreshold)){
    #overlay a mincount filter
    idxmin <- counts >= countThreshold
    frac <- rowSums(idxmin)/ncol(idxmin)
    idx <- frac >= sampleFraction
    counts <- counts[idx,]
    if (xClass == "DGEobj") x <- x[idx,]
    genelength <- genelength[idx]

    if (verbose == TRUE)
      tsmsg(stringr::str_c(sum(idx), " of ", length(idx), " genes retained by the low count filter."))
  }

  if (xClass == "DGEobj") return(x) else return(counts)

}



