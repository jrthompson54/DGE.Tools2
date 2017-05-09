### Function isoformFrac ###
#' Function isoformFrac
#'
#' Takes a DGEobj as input  and adds an item containing isoform fraction data.
#'
#' Isoform Fraction is calculated using TPM data, as length normalized data
#' is required since different isoforms have different total exon lengths.  Isoform
#' fraction is calculated simply as the isoform TPM divided by the total TPM
#' for all isoforms of a given gene.
#' 
#' TPM is calculated from counts using the ExonLength column from the gene/isoform
#' data.  If normalization is specified, TMM or other selected normalization
#' is applied to CPM data (using edgeR functions) before converting to TMM.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, DGEobj
#'
#' @param dgeObj  A isoform level dgeObj created by function initDGEobj
#' @param normalize Default = TRUE and invokes TMM normalization. FALSE uses
#' TPM without TMM normalization.  Other allowed values are: "TMM","RLE","upperquartile". 
#' Invokes edgeR::calcNormFactors for normalization.
#'
#' @return A DGEobj with an isoform fraction dataframe added
#'
#' @examples
#'    MyDgeObj <- isoformFrac(MyDgeObj)
#'
#' @import dplyr tidyr DGEobj magrittr assertthat
#'
#' @export
isoformFrac <- function(dgeObj){

    assertthat::assert_that(class(dgeObj)[[1]] == "DGEobj",
                          attr(dgeObj, "level") == "isoform")

    #calculate sum of isoform TPM for each gene and sample
    counts <- getItem(dgeObj, "counts")
    tData <- getItem(dgeObj, "isoformData")
    tpm <- convertCounts(counts, unit="TPM", geneLength=tData$ExonLength) %>%
      as.data.frame
    tpm$GeneID <- tData$GeneID
    tpm$TranscriptID <- rownames(tpm)

    #Calculate isoform fraction
    tpmtidy <- tidyr::gather(tpm, key=sample, value=tpm, -GeneID, -TranscriptID) %>%
      dplyr::group_by (sample, GeneID) %>%
      dplyr::mutate (geneTPM = sum(tpm),
              isofrac = tpm/geneTPM)

    #drop uneeded columns
    tpmtidy$tpm <- NULL
    tpmtidy$geneTPM <- NULL

    #now spread to an isoformPct matrix
    IsoformFrac <- spread(tpmtidy, sample, isofrac) %>% as.data.frame
    #set rownames to transcript ID and remove ID columns
    rownames(IsoformFrac) <- IsoformFrac$TranscriptID
    IsoformFrac$GeneID <- NULL
    IsoformFrac$TranscriptID <- NULL
    
    #debug
    # saveRDS(IsoformFrac, "isoformFraction.RDS")
    
    #add isoform fraction to assays.
    funArgs <- match.call()
    dgeObj <- addItem(dgeObj, IsoformFrac, itemName="isoformFrac",
                      itemType="assay", funArgs=funArgs,
                      parent="counts")

    
    return(dgeObj)
}

