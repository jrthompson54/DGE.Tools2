### Function convertCounts ###
#' Function convertCounts
#'
#' Takes a count matrix as input and converts to other desired units.  Supported
#' units include CPM, FPKM, FPK and TPM.  Output units can be logged
#' and/or normalized.  Calculations are performed using edgeR functions except
#' for the conversion to TPM which is converted from FPKM using the formula provided
#' by [Harold Pimental](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
#' 
#' geneLength is a vector where length(geneLength) == nrow(counts).  If you pass
#' an RSEM effectiveLength matrix, rowMeans(effectiveLength) is used.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, unit conversion
#'
#' @param counts A numeric matrix or dataframe of N genes x M Samples.  All columns
#' must be numeric.
#' @param unit  Required. One of CPM, FPKM, FPK or TPM.
#' @param geneLength Required for length-normalizes units (TPM, FPKM or FPK)
#' @param log Default = FALSE.  Set TRUE to return Log2 values. Log conversion
#' employs the edgeR method which uses an average prior of 0.25 moderated by the
#'    library size.
#' @param normalize Default = FALSE. TRUE activates TMM normalization.
#'    Other allowed values are: "TMM","RLE","upperquartile".  Invokes
#'    edgeR::calcNormFactors for normalization.
#'
#' @return A matrix in the new unit space
#'
#' @examples
#' #TMM normalized Log2TPM
#' Log2TPM = convertCounts(mycounts),
#'                       unit="TPM",
#'                       geneLength=gene.annotation$ExonLength,
#'                       log=TRUE,
#'                       normalize=TRUE)
#'
#' #un-normalized CPM (not logged)
#' RawCPM = convertCounts(MyCounts,
#'                       unit="CPM",
#'                       log=FALSE,
#'                       normalize=FALSE)
#'
#' @import edgeR zFPKM dplyr magrittr assertthat
#'
#' @export
convertCounts <- function(counts,
                          unit,
                          geneLength,
                          log=FALSE,
                          normalize=FALSE,
                          PlotDir,
                          PlotFile,
                          FacetTitles) {
#   Good explanation of prior values in edgeR vs. voom CPM/elist values:
#   https://support.bioconductor.org/p/59846/

    .calcCPM <- function(counts, log, normalize){
        if (nrow(counts) < 10000)
            warning('You should use the whole dataset when calculating CPM, not a subset.')
        counts %>%
        DGEList %>%
        calcNormFactors(method=normalize) %>%
        cpm(log=log)
    }

    .calcFPKM <- function(counts, log, normalize, geneLength){
        if (nrow(counts) < 10000)
            warning('You should use the whole dataset when calculating FPKM, not a subset.')
        counts %>%
        DGEList %>%
        calcNormFactors(method=normalize) %>%
        rpkm(log=log, gene.length=geneLength)
    }

    .calcTPM <- function(counts, log, normalize, geneLength){
        if (nrow(counts) < 10000)
            warning('You should use the whole dataset when calculating TPM, not a subset.')
        # fpkm <- counts %>%
        #     DGEList %>%
        #     calcNormFactors(method=normalize) %>%
        #     rpkm(log=log, gene.length=geneLength)
        fpkm <- .calcFPKM(counts, log=FALSE, normalize=normalize, geneLength = geneLength)
        if (log==FALSE)
            TPM <- .fpkmToTpm(fpkm)
        else
            TPM <- log2(.fpkmToTpm(2^fpkm))
        return(TPM)
    }

    .calcFPK <- function(counts, log, normalize, geneLength){
        if (tolower(normalize) == 'none'){
            #check for zero geneLength just in case
            if (min(geneLength) == 0) geneLength <- geneLength + 1

            FPK <- counts / (geneLength/1000)

            if (log == TRUE)
                FPK <- log2(FPK + 0.5)
        } else {
            #why would you want normalized FPK???
            #See Mark Robinson response on how to get normalized counts
            #from edgeR
            #https://stat.ethz.ch/pipermail/bioconductor/2012-July/046795.html
            #JRT: I haven't validated this because I can think of no reason
            #why you'd want normalized FPK
            stop("FPK should not be normalized")
        }
        return(FPK)
    }

    .fpkmToTpm <- function(fpkm, warn=TRUE)
    {
      exp(log(fpkm) - log(colSums(fpkm, na.rm=TRUE)) + log(1e6))
    }

### MAIN ###
#counts and unit are absolutely required
if (missing(counts))
    stop ("counts is a required argument")
    
if (missing(unit))
    stop ("unit is a required argument")
    
if (toupper(unit) %in% c('FPKM', 'TPM', 'FPK')){ 
    #in these cases geneLength is required
    if (missing(geneLength))
        stop("geneLength is required for unit = FPK|FPKM|TPM")
    
    if (class(geneLength)[[1]] == "matrix") #flatten to a vector
        geneLength <- rowMeans(geneLength, na.rm=TRUE)
}

#Coerce counts to a matrix
result <- try({counts <- as.matrix(counts)}, silent=TRUE)
if (class(result)[[1]] == "try-error")
    stop("Couldn't coerce counts to a numeric matrix!")

#Make sure geneLength is correct length
if (!missing(geneLength))
    if (length(geneLength) != nrow(counts))
        stop('Length(geneLength) does not match rowcount of counts')

#set defaults
if (missing(log))
    log=FALSE
if (missing(normalize))
    normalize='none'
if (is.logical(normalize)) #convert to char
    if (normalize == TRUE)
        normalize <- 'TMM'
    else normalize <- 'none'

result <- switch(toupper(unit),
       "CPM" = .calcCPM(counts, log, normalize),
       "FPKM" = .calcFPKM(counts, log, normalize, geneLength),
       "FPK" = .calcFPK(counts, log, normalize, geneLength),
       "TPM" = .calcTPM(counts, log, normalize, geneLength)
)
return(result)
}

### Function tpm ###
#' Function tpm
#'
#' Takes a DGEObj as input and converts counts_orig to tpm units.  
#' For GeneLength in the calculation, it takes geneData$ExonLength (for Omicsoft data) 
#' or uses rowMeans(effectiveLength) for data derived from Xpress.
#' 
#' Internally it used edgeR::fpkm to calculate fpkm and converts to tpm
#'  using the formula provided  by [Harold Pimental](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, unit conversion
#'
#' @param counts A numeric matrix or dataframe of N genes x M Samples.  All columns
#' must be numeric.
#' @param unit  Required. One of CPM, FPKM, zFPKM, FPK or TPM.
#' @param geneLength Required for length-normalizes units (TPM, FPKM, zFPKM or FPK)
#' @param log Default = FALSE.  Set TRUE to return Log2 values. Log conversion
#' employs the edgeR method which uses an average prior of 0.25 moderated by the
#'    library size.
#' @param normalize Default = FALSE. TRUE activates TMM normalization.
#'    Other allowed values are: "TMM","RLE","upperquartile".  Invokes
#'    edgeR::calcNormFactors for normalization.
#' @param PlotDir Only applies to zFPKM. Specifies folder for zFPKM distribution plot.
#' @param PlotFile Only applies to zFPKM. Default = "zFPKM_Distribution.png"
#' @param FacetTitles Only applies to zFPKM. Turn facet titles on or off
#'    (Default = TRUE)
#'
#' @return A matrix in the new unit space
#'
#' @examples
#' myTPM <- tpm(dgeobj) 
#'
#' @import DGEobj dplyr magrittr assertthat
#'
#' @export
tpm <- function(dgeo){
    
    assertthat::assert_that(class(dgeo)[[1]] == "DGEobj")

    #default to gene level
    level <- "gene"
    if (attr(dgeo, "level") == "isoform")
        level <- "isoform"
    
    #Not supporting Exon level yet.  Not sure where to grab the geneLength data from for exon level
    if (attr(dgeo, "level") == "exon")
        stop("Exon level not supported; contact JRT and make the request")
    
    if (level == "gene") 
        rowdata <- getItem(dgeo, "geneData_orig")
    else 
        rowdata <- getItem(dgeo, "isoformData_orig")
    
    #get genelength depending on source data (omicsoft or xpress)
    if (attr(dgeo, "source") == "Omicsoft"){ #Omicsoft data
        geneLength <- rowdata$ExonLength
    } else if ("effectiveLength_orig" %in% names(dgeo)) { #use rowMeans(effectiveLength)
        geneLength <- rowMeans(getItem(dgeo, "effectiveLength_orig"), na.rm=TRUE)
    } else {
        stop("Couldn't find gene/isoform length data in DGEobj")
    }

    TPM <- convertCounts(getItem(dgeo, "counts_orig"), geneLength=geneLength, unit="tpm",
                         log=FALSE, normalize=FALSE)
        
}
