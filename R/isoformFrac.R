### Function isoformFrac ###
#' Function isoformFrac
#'
#' Takes a DGEobj as input (transcript level data) and adds an assay item containing isoform fraction data.
#'
#' Isoform Fraction is calculated using length normalized data (FPKM or TPM), as
#' length normalized data is required because different isoforms have different
#' total exon lengths.  If FPKM is specified, you can also specify a
#' normalization (via edgeR::calcNormFactors). Isoform fraction is calculated
#' simply as the isoform intensity divided by the total gene intensity for all isoforms of
#' a given gene.
#'
#' TPM or FPKM are calculated directly from counts using all data in the dgeObj.
#' I recommend performing low intensity gene filtering before running isoformFrac.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, DGEobj
#'
#' @param dgeObj  An isoform level DGEobj created by function initDGEobj,
#'   Xpress2DGEO or OmicsoftToDgeObj.  Counts and isoformData must be present
#'   in the DGEobj (required).  isoformData$ExonLength must be present or assay =
#'   "effectiveLength" must be present.
#' @param dataType One of "fpkm" or "tpm" (default="fpkm")
#' @param normalize Default = "TMM" and invokes TMM normalization. Other allowed
#'   values are: "RLE","upperquartile", "none". Invokes edgeR::calcNormFactors for
#'   normalization.  Only invoked when dataType="fpkm".  This is because
#'   applying TPM essentially erases any prior column scaling so TMM and similar
#'   normalizations have no effect.
#'
#' @return A DGEobj with an isoform fraction dataframe added
#'
#' @examples
#'    MyDgeObj <- isoformFrac(MyDgeObj)
#'
#' @import magrittr
#' @importFrom dplyr group_by mutate
#' @importFrom assertthat assert_that
#' @importFrom tidyr gather spread
#' @importFrom DGEobj getItem addItem
#'
#' @export
isoformFrac <- function(dgeObj, dataType="fpkm", normalize="tmm"){

    assertthat::assert_that(class(dgeObj)[[1]] == "DGEobj",
                          attr(dgeObj, "level") == "isoform")

    #calculate sum of isoforms for each gene and sample
    counts <- DGEobj::getItem(dgeObj, "counts")
    isoformData <- DGEobj::getItem(dgeObj, "isoformData")

    #Xpress support:
    #If present, Need to take rowMeans of assay=effectiveLength and create
    #isoformData$ExonLength.
    if (tolower(attr(dgeObj, "source")) == "xpress"){#It's an Xpress dataset
      efflen <- DGEobj::getItem(dgeObj, "effectiveLength")
      if (!is.null(efflen)){
        isoformData <- DGEobj::getItem(dgeObj, "isoformData")
        isoformData$ExonLength <- rowMeans(efflen, na.rm=TRUE)
      } else {
        stop ("Effective Length data not found in Xpress DGEobj!")
      }
    }

    omicData <- switch (tolower(dataType),
      "fpkm" = convertCounts(counts, unit="fpkm",
                             geneLength=isoformData$ExonLength,
                             normalize = normalize),
      "tpm" = convertCounts(counts, unit="TPM",
                            geneLength=isoformData$ExonLength,
                            normalize = "none")
    ) %>% as.data.frame()

    omicData$GeneID <- isoformData$GeneID
    omicData$TranscriptID <- rownames(omicData)

    #Calculate isoform fraction
    omictidy <- tidyr::gather(omicData, key=sample, value=intensity, -GeneID, -TranscriptID) %>%
      dplyr::group_by (sample, GeneID) %>%
      dplyr::mutate (geneTotal = sum(intensity),
              isofrac = intensity/geneTotal)

    #drop uneeded columns
    omictidy$intensity <- NULL
    omictidy$geneTotal <- NULL

    #now spread to an isoformPct matrix
    IsoformFrac <- spread(omictidy, sample, isofrac) %>% as.data.frame
    #set rownames to transcript ID and remove ID columns
    rownames(IsoformFrac) <- IsoformFrac$TranscriptID
    IsoformFrac$GeneID <- NULL
    IsoformFrac$TranscriptID <- NULL

    #put isoform fraction data in same order as dgeobj
    IsoformFrac <- IsoformFrac[rownames(dgeObj),]

    #debug
    # saveRDS(IsoformFrac, "isoformFraction.RDS")

    #add isoform fraction to assays.
    funArgs <- match.call()
    dgeObj <- DGEobj::addItem(dgeObj, IsoformFrac,
                      itemName=paste("isoformFrac", dataType, sep="_"),
                      itemType="assay", funArgs=funArgs,
                      parent="counts")

    return(dgeObj)
}

