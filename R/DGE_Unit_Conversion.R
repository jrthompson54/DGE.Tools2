#DGE Unit conversion functions
# from: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

### Function countToTpm ###
#' Function countToTpm:  Convert Counts to TPM
#'
#' To get a proper TPM value you should use the unfiltered counts matrix as input to this function.
#' This is because a column sum is involved in the TPM calculation and the column sum should include all mapped reads.
#' 
#' Reference = https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#' Modified by JRT to work with an input matrix instead of just a vector.
#'
#' @author https://haroldpimentel.wordpress.com/
#' @keywords read counts; TPM
#'
#' @param counts a matrix of genes(rows) and samples(col)x = assay
#' @param effLen vector or gene lengths (length must equal length(counts)).
#'
#' @return matrix of TPM values
#'
#' @examples
#' MyTPM = countToTpm(counts, MyCounts)
#'
#' @export
countToTpm <- function(counts, effLen)
{

	 if (nrow(counts) < 20000) {
		print ("WARNING: The TPM calculation will not be accurate unless you use the full unfiltered count matrix")
	 }
    rate <- log(counts) - log(effLen)
    denom <- log(colSums(exp(rate)))
    exp(rate - denom + log(1e6))
    
}

### Function countToFpkm ###
#' Function countToFpkm:  Convert Counts matrix to TPM
#'
#' To get a proper FPKM value you should use the unfiltered FPKM matrix as input to this function.
#' This is because a column sum is involved in the TPM calculation and the column sum should include all mapped reads.
#' 
#' Reference = https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#' Modified by JRT to work with an input matrix instead of just a vector.
#'
#' @author https://haroldpimentel.wordpress.com/
#' @keywords read counts; TPM
#'
#' @param counts a matrix of genes(rows) and samples(col)
#' @param effLen vector or gene lengths (length must equal length(counts)).
#'
#' @return matrix of FPKM values
#'
#' @examples
#' MyFPKM = countToTpm(counts, MyCounts)
#'
#' @export
countToFpkm <- function(counts, effLen)
{
	 if (nrow(counts) < 20000) {
		print ("WARNING: The FPKM calculation will not be accurate unless you use the full unfiltered FPKM matrix")
	 }
    N <- colSums(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

### Function fpkmToTpm ###
#' Function fpkmToTpm:  Convert FPKM to TPM
#'
#' To get a proper TPM value you should use the unfiltered FPKM matrix as input to this function.
#' This is because a column sum is involved in the TPM calculation and the column sum should include all mapped reads.
#' 
#' Reference = https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#' Modified by JRT to work with an input matrix instead of just a vector.
#'
#' @author https://haroldpimentel.wordpress.com/
#' @keywords read counts; TPM; FPKM
#'
#' @param fpkm a matrix of genes(rows) by samples(col)
#'
#' @return matrix of TPM values
#'
#' @examples
#' MyTPM = fpkmToTpm(MyFPKM)
#'
#' @export
fpkmToTpm <- function(fpkm, warn=TRUE)
{
	 if (nrow(fpkm) < 20000 & warn==TRUE) {
		print ("WARNING: The TPM calculation will not be accurate unless you use the full unfiltered fpkm matrix")
	 }
    exp(log(fpkm) - log(colSums(fpkm)) + log(1e6))
}

### Function Log2CPM_To_Log2FPKM ###
#' Function Log2CPM_To_Log2FPKM:  Convert Log2CPM to Log2FPKM
#'
#' Convert Log2CPM from a voom elist (elist$E) to LOG2FPKM
#' to provide a Log2FPKM value to compare between
#' genes that is fully normalized and voom adjusted.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords CPM; FPKM
#'
#' @param log2cpm A dataframe or matrix of Log2CPM values typically
#' from the voom  Elist ($E)
#' @param exonLength vector of total exon length from the design
#' (rowRanges) data.
#'
#' @return matrix of FPKM values
#'
#' @examples
#' MyFPKM = cpmToFpkm((2^log2cpm), exonlength)
#'
#' @export
Log2CPM_To_Log2FPKM <- function (log2cpm, exonLength) {
  #FPKM = CPM / (exonLength/1000)
  log2fpkm = log2cpm - log2(exonLength/1000)
}

### Function cpmToFpkm ###
#' Function cpmToFpkm:  Convert Counts to TPM
#'
#' Reference = https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#'
#' @author https://haroldpimentel.wordpress.com/
#' @keywords CPM; FPKM
#'
#' @param counts a matrix of genes(rows) and samples(col)
#' @param effLen vector or gene lengths (length must equal length(counts)).
#'
#' @return matrix of FPKM values
#'
#' @examples
#' MyFPKM = cpmToFpkm(MyCPM, exonlengths)
#'
#' @export
cpmToFpkm <- function (CPM, exonLength) {
  FPKM = CPM / (exonLength/1000)
}

### Function counts2TMM_TPM ###
#' Function counts2TMM_TPM:  Convert Counts to TMM normalized TPM
#'
#' TMM alone does not capture the compositional bias normalization provided by TMM normalization.
#' This function uses edgeR functions to calculate TMM normalized CPM and then coverts it to TMM units.
#' Note that you should use a complete count matrix of all genes for the TPM units to be correct as
#' the conversion to TPM involves a column sum.
#' 
#' Reference = https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#' Modified by JRT to work with an input matrix instead of just a vector.
#'
#' @author https://haroldpimentel.wordpress.com/
#' @keywords read counts; TPM; FPKM
#'
#' @param counts a matrix of genes(rows) by samples(col); use all genes, not a filtered gene list
#' @param genelength a required vector of gene lengths in the same order as rows of counts.  This is typically the effective length
#' 	from the RSEM calculation.  You can get this from the "exonlength" column in Omicsoft Annotation and the effective
#'    length data should also be available for Xpress projects as of Aug2016.
#'
#' @return matrix of TMM normalized TPM values
#'
#' @examples
#' MyTPM = counts2TMM_TPM(MyCounts, MyGeneLengths)
#'
#' @import magrittr edgeR
#'
#' @export
counts2TMM_TPM <- function(counts, genelength)
{
	 if (nrow(counts) < 20000) {
		print ("WARNING: The TPM calculation will not be accurate unless you use the full unfiltered count matrix")
	 }

	TPM <- counts %>% DGEList %>% calcNormFactors(method = "TMM") %>%
	                rpkm(log=FALSE, gene.length=genelength) %>% fpkmToTpm(warn=FALSE) 
}

