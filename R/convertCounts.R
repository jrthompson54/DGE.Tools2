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
#' an RSEM effectiveLength matrix, rowMeans(effectiveLength) is used (because edgeR functions
#' only accept a vector for effectiveLength).
#'
#' Note that log2 values for CPM, TPM and FPKM employ edgeR's prior.count handling to avoid divide by zero.
#'
#' @author John Thompson, \email{jrt@@thompsonclan.org}
#' @keywords RNA-Seq, unit conversion
#'
#' @param counts A numeric matrix or dataframe of N genes x M Samples.  All columns
#' must be numeric.
#' @param unit  Required. One of CPM, FPKM, FPK or TPM.
#' @param geneLength A vector or matrix of gene lengths. Required for length-normalizes units (TPM, FPKM or FPK).
#'    If you supply a matrix, the rowMeans are calculated and used.
#' @param log Default = FALSE.  Set TRUE to return Log2 values.
#'    Employs edgeR functions which use an prior.count of 0.25 scaled by the library size.
#' @param normalize Default = "none". Other options: "TMM", "RLE", "upperquartile"
#'  Invokes edgeR::calcNormFactors for normalization. Upperquartile uses the 75th percentile.  Normalize settings are case insensitive.
#' @param prior.count Average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE. (Default dependent on method;
#'   0 for TPM, 0.25 for CPM and FPKM)
#'   The prior.count is passed to edgeR cpm and rpkm functions and applies to logTPM, logCPM and logFPKM calculations.
#'
#' @return A matrix in the new unit space
#'
#' @examples
#' #TMM normalized Log2FPKM
#' Log2FPKM = convertCounts(mycounts),
#'                       unit="fpkm",
#'                       geneLength=gene.annotation$ExonLength,
#'                       log=TRUE,
#'                       normalize="tmm")
#'
#' #un-normalized CPM (not logged)
#' RawCPM = convertCounts(MyCounts,
#'                       unit="CPM",
#'                       log=FALSE,
#'                       normalize="none")
#'
#' @import magrittr
#' @importFrom edgeR cpm rpkm expandAsMatrix calcNormFactors DGEList
#' @importFrom assertthat assert_that
#' @importFrom glue glue
#'
#' @export
convertCounts <- function(counts,
                          unit,
                          geneLength,
                          log=FALSE,
                          normalize="none",
                          prior.count=NULL,
                          debug=FALSE) {
    #   Good explanation of prior values in edgeR vs. voom CPM/elist values:
    #   https://support.bioconductor.org/p/59846/

    ### MAIN ###
    #counts and unit are absolutely required
    if (missing(counts))
        stop ("counts is a required argument")

    if (missing(unit))
        stop ("unit is a required argument")

    unit <- toupper(unit)
    if (unit %in% c('FPKM', 'TPM', 'FPK')){
        #in these cases geneLength is required
        if (missing(geneLength))
            stop("geneLength is required for unit = FPK|FPKM|TPM")

        if (class(geneLength)[[1]] == "matrix") #flatten to a vector
            geneLength <- rowMeans(geneLength, na.rm=TRUE)
    }
    #make normalize method case insensitive (calcNormFactors is case sensitive)
    if (toupper(normalize) %in% c("TMM", "RLE")) #these have to be uppercase
        normalize <- toupper(normalize)
    if (toupper(normalize) %in% c("UPPERQUARTILE", "NONE")) #these have to be lowercase
        normalize <- tolower(normalize)
    if (toupper(normalize) %in% c("TMMWZP"))
      normalize <- "TMMwzp"

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
    if (is.logical(normalize)){ #don't encourage logicals; here for backward compatibility
        if (normalize == TRUE)
            normalize <- 'TMM'
        if (normalize == FALSE)
            normalize <- 'none'
    }

    if (is.null(prior.count)) {
        if (log == FALSE){
            prior.count <- 0 #not used when log=F
        }
        else if (unit == "TPM") {
            prior.count <- 0
        }
        else {
            prior.count <- 0.25
        }
    }

    result <- switch(toupper(unit),
                     "CPM" = calcCPM(counts, log, normalize, prior.count, debug),
                     "FPKM" = calcFPKM(counts, log, normalize, geneLength, prior.count, debug),
                     "FPK" = calcFPK(counts, log, normalize, geneLength, prior.count, debug),
                     "TPM" = calcTPM(counts, log, normalize, geneLength, prior.count, debug)
    )
    return(result)
}

############### Helper Functions

calcCPM <- function(counts, log, normalize, prior.count, debug){
  #debug
  if (debug) print(glue("CPM; log= {log}; normalize={normalize}; prior={prior.count}"))
  if (nrow(counts) < 10000)
    warning('You should use the whole dataset when calculating CPM, not a subset.')
  counts %>%
    edgeR::DGEList() %>%
    edgeR::calcNormFactors(method=normalize) %>%
    edgeR::cpm(log=log, prior.count=prior.count)
}

calcFPKM <- function(counts, log, normalize, geneLength, prior.count, debug){
  #debug
  if (debug) print(glue("FPKM; log= {log}; normalize={normalize}; prior={prior.count}"))
  if (nrow(counts) < 10000)
    warning('You should use the whole dataset when calculating FPKM, not a subset.')
  counts %>%
    edgeR::DGEList() %>%
    edgeR::calcNormFactors(method=normalize) %>%
    edgeR::rpkm(log=log, gene.length=geneLength, prior.count=prior.count)
}

calcTPM <- function(counts, log, normalize, geneLength, prior.count, debug){
    #debug
    if (debug) print(glue("TPM; log= {log}; normalize={normalize}; prior={prior.count}"))
    if (nrow(counts) < 10000)
        warning('You should use the whole dataset when calculating TPM, not a subset.')
    if (normalize != "none")
        warning(glue('TPM normalization overides {normalize} normalization!'))
    if (prior.count != 0 && log == TRUE)
        warning("Using a prior.count for logtpm calculations is not recommended and may produce unpredictable results!")
    fpkm <- calcFPKM(counts, log=log, normalize=normalize,
                     geneLength = geneLength, prior.count=prior.count, debug)

    if (log==FALSE){
        TPM <- fpkmToTpm(fpkm)
    } else {
        TPM <- log2(fpkmToTpm(2^fpkm))
    }
    return(TPM)
}

calcFPK <- function(counts, log, normalize, geneLength, prior.count, debug){
    #debug
    if (debug) print(glue("FPK; log= {log}; normalize={normalize}; prior={prior.count}"))
    if (tolower(normalize) == 'none'){
        #check for zero geneLength just in case
        if (min(geneLength) == 0) geneLength <- geneLength + 1

        FPK <- counts / (geneLength/1000)

        if (log == TRUE)
            FPK <- log2(FPK + prior.count)
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

fpkmToTpm <- function(fpkm)
{
    colSumMat <- edgeR::expandAsMatrix(colSums(fpkm, na.rm=TRUE), byrow=TRUE, dim=dim(fpkm))
    #exp(log(fpkm) - log(colSumMat) + log(1e6)) #only works with fpkm vector
    fpkm / colSumMat * 1e6
}

### Function tpm.on.subset ###
#' Function tpm.on.subset
#'
#' Use this to calculate TPM for a heavily subsetted DGEobj.  It calculates TPM
#' using the original data but returns a DGEobj with the subset.
#'
#' Takes a DGEObj as input. Uses all data (pre-genefiltering; counts_orig) to
#' calculate TPM values. Then optionally filters any genes that were already
#' filtered out.
#'
#' For GeneLength in the calculation, it takes geneData$ExonLength (for Omicsoft data)
#' or uses rowMeans(effectiveLength) for data derived from Xpress.
#'
#' Internally it uses edgeR::fpkm to calculate fpkm and converts to tpm
#' using the formula provided  by [Harold Pimental](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
#'
#' @author John Thompson, \email{jrt@@thompsonclan.org}
#' @keywords RNA-Seq, unit conversion
#'
#' @param dgeo A DGEobj data structure
#' @param applyFilter If TRUE, reduces to the filtered gene list. FALSE returns
#'   all genes in the raw data. [Default = TRUE]
#'
#' @return A matrix in the new unit space
#'
#' @examples
#' myTPM <- tpm(dgeobj)
#'
#' @import DGEobj magrittr
#' @importFrom assertthat assert_that
#' @export
tpm.on.subset <- function(dgeo, applyFilter=TRUE){

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

    #remove filtered out genes
    if (applyFilter == TRUE){
        idx <- rownames(TPM) %in% rownames(getItem(dgeo, "counts"))
        TPM <- TPM[idx,]
    }
    return(TPM)
}

### Function tpm.direct ###
#' Function tpm.direct
#'
#' Takes a counts and genelength as input and converts to tpm units using the equation from
#' [Harold Pimental](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
#'
#' The result should be the same as using convertCounts with normalize='tpm' and log=FALSE
#'
#' Genelength can be a vector (length == nrow(counts) or a matrix (same dim as counts).
#' The genelength is used as is, or optionally collapsed to a vector by rowMeans.
#'
#' @author John Thompson, \email{jrt@@thompsonclan.org}
#' @keywords RNA-Seq, unit conversion
#'
#' @param counts A numeric matrix of N genes x M Samples.  All columns
#' must be numeric.
#' @param geneLength A vector or matrix of gene lengths.
#' @param collapse T/F determines whether to use rowMeans on the genelength matrix [Default = FALSE]
#'
#' @return A matrix of TPM values
#'
#' @examples
#' myTPM <- tpm.classic(mycounts, mygenelength)
#'
#' @import magrittr
#' @importFrom edgeR expandAsMatrix
#' @importFrom assertthat assert_that
#' @export
tpm.direct <- function (counts, genelength, collapse=FALSE){
    #equation for TPM from https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

    if (!is.matrix(counts)){ #try to coerce
        result <- try(counts <- as.matrix(counts), silent=TRUE)
        if (class(result) == "try-error")
            stop("Couldn't coerce counts to a matrix!")
    }

    if (is.vector(genelength)){
        assertthat::assert_that(length(genelength) == nrow(counts))
    } else { #try to coerce
        if (!is.matrix(genelength)){
            result <- try(genelength <- as.matrix(genelength), silent=TRUE)
            if (class(result) == "try-error")
                stop("Couldn't coerce genelength to a matrix!")
            assertthat::assert_that(all(dim(counts) == dim(genelength)))
        }
    }

    if (collapse & is.matrix(genelength)) { #reduce genelength to a vector
        genelength <- rowMeans(genelength, na.rm=TRUE)
    }

    #the calculation  (fpk / colsum(fpk) ) * 10e6
    fpb <- counts / genelength
    sumfpb <- colSums(fpb)
    tpm <- fpb / edgeR::expandAsMatrix(sumfpb, byrow = TRUE, dim=dim(fpb)) * 1e6
    # tpm <- t(t(fpb)/sumfpb) * 1e6
}
