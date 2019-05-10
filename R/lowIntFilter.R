### Function lowIntFilter ###
#' Function  lowIntFilter
#'
#' Takes a DGEobj as input and applies a combination of low
#' intensity filters. raw count, zFPKM, TPM and/or FPK filters are supported.  A gene
#' must pass all active filters.  Not setting a threshold argument inactivates that threshold.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq; counts; low intensity
#'
#' @param x A DGEobj with RNA-Seq (counts) data
#'   (required).
#' @param countThreshold Genes below this threshold are removed (recommend 10).
#'   Set to 0 to disable this filter.
#' @param zfpkmThreshold  Genes below this threshold are removed (-3.0 is recommended)
#' @param fpkThreshold Genes below this threshold are removed (recommend 5).
#' @param tpmThreshold Genes below this threshold are removed (tpm is supported
#'   by request, but FPK is a better length-normalized value to use as a filter)
#' @param sampleFraction The proportion of samples that must meet the thresholds
#'   (Default = 0.5). Range >0 and <=1.
#' @param genelength Vector of genelengths for rows of x. Required for FPK and
#'   zFPKM filters, unless x is a DGEobj.  If a DGEobj is supplied, genelength is
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
#' @importFrom DGEobj getItem
#'
#' @export
lowIntFilter <- function(x, zfpkmThreshold, fpkThreshold, countThreshold, tpmThreshold,
                         sampleFraction=0.5, genelength=NULL,
                         verbose=FALSE)
{

  assertthat::assert_that("DGEobj" %in% class(x))

  #can't use both tpm and zfpkm
  if (!missing(zfpkmThreshold) & !missing(tpmThreshold)){
    stop("Must use zfpkmThreshold or tpmThreshold but not both!")
  }

  #need counts
  counts <-getItem(x, "counts")
  genelength <- x$geneData$ExonLength

  starting_rowcount <- nrow(counts)

  #get genelength from dgeobj
  if (is.null(genelength))  #user supplied genelength supercedes genelength from DGEobj
    genelength <-x$geneData$ExonLength

  #apply zFPKM threshold
  if (!missing(zfpkmThreshold)){
    assertthat::assert_that(!is.null(genelength),
                            length(genelength) == nrow(x$counts))
    fpkm <- convertCounts(x$counts, unit="fpkm", geneLength = genelength) #class(fpkm) = "matrix"

    #Xpress counts produce rows filled with NaNs.  Need to filter these out first before calculating zFPKM.
    idx <- complete.cases(fpkm)
    x <- x[idx,]
    fpkm <- fpkm[idx,]
    genelength <- genelength[idx]
    nan_genes <- sum(!idx)
    if (nan_genes > 0 & verbose==TRUE)
      message(stringr::str_c(nan_genes, " genes with NaN FPKM values removed."))

    #calculate zfpkm
    zfpkm <-as.matrix(zFPKM::zFPKM(as.data.frame(fpkm)))

    #create index for zfpkm >= zfpkmThreshold in fracThreshold of samples
    idx_zfpkm <- zfpkm >= zfpkmThreshold
    frac <- rowSums(idx_zfpkm) / ncol(idx_zfpkm)
    fpkmidx <- frac >= sampleFraction
    x <- x[fpkmidx,]
    genelength <- genelength[fpkmidx]

    if (verbose == TRUE) {
      message(stringr::str_c(sum(fpkmidx), " of ", starting_rowcount, " genes retained by the zFPKM filter."))
    }
  }

  #apply TPM threshold
  if (!missing(tpmThreshold)){
    assertthat::assert_that(!is.null(genelength),
                            length(genelength) == nrow(counts))
    tpm <- convertCounts(x$counts, unit="tpm", geneLength = genelength)

    #Xpress counts produce rows filled with NaNs.  Need to filter these out first before calculating tpm
    idx <- complete.cases(tpm)
    x <- x[idx,]
    tpm <- tpm[idx,]
    genelength <- genelength[idx]
    nan_genes <- sum(!idx)
    if (nan_genes > 0 & verbose==TRUE)
      message(stringr::str_c(nan_genes, " genes with NaN TPM values removed."))

    #create index for tpm >=tpmThreshold in fracThreshold of samples
    idx_tpm <- tpm >= tpmThreshold
    frac <- rowSums(idx_tpm) / ncol(idx_tpm)
    tpmidx <- frac >= sampleFraction
    x <- x[tpmidx,]
    genelength <- genelength[tpmidx]

    if (verbose == TRUE) {
      message(stringr::str_c(sum(tpmidx), " of ", starting_rowcount, " genes retained by the TPM filter."))
    }
  }


  #apply FPK threshold
  if (!missing(fpkThreshold)){
    assertthat::assert_that(!is.null(genelength),
                            length(genelength) == nrow(x))
    fpk <- convertCounts(x$counts, unit="fpk", geneLength = genelength)
    #keep FPK >= fpkThreshold in fracThreshold of samples
    idxfpk <- fpk >= fpkThreshold
    frac <- rowSums(idxfpk)/ncol(idxfpk)
    idx <- frac >= sampleFraction
    x <- x[idx,]
    genelength <- genelength[idx]

    if (verbose == TRUE)
      message(stringr::str_c(sum(idx), " of ", length(idx), " genes retained by the FPK filter."))
  }

  #apply count threshold
  if (!missing(countThreshold)){
    #overlay a mincount filter
    idxmin <- x$counts >= countThreshold
    frac <- rowSums(idxmin)/ncol(idxmin)
    idx <- frac >= sampleFraction
    x <- x[idx,]
    genelength <- genelength[idx]

    if (verbose == TRUE)
      message(stringr::str_c(sum(idx), " of ", length(idx), " genes retained by the low count filter."))
  }

  return(x)

}



