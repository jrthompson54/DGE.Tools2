# DGE_Tools.R
# Author: JRT 28Dec2015
#
# A set of functions to facilitate differential gene expression analysis
# of data exported from the Omicsoft RNA-Seq Pipeline.  Easily Adaptable
# to other data sources of, minimally, count data, row(gene) annotation
# and col(sample) annotation.
#
# Main Functions:
# Build_RSE:  Read in a named list of expresson data files
# plotDispersion: generate an edgeR dispersion plot
# DGE_Calc: run edgeR/voom pipeline and run a fit
#
# Utility Functions: not exported
# Txt2DF: Utility function called by Build_EO
# fpkmToTpm: Utility function called by Build_RSE
# tsmsg: a timestamped message
# df2gr: Convert a DF with chr, start, end, strand column to a GRanges object.

###LIBRARIES used
# library(magrittr)
# library(gdata)
# library(openxlsx)
# library(GenomicRanges)
# library(SummarizedExperiment)
# library(dplyr)
# library(reshape2)
# library(ggplot2)

### DATA STRATEGY ###
# DGE analysis will involve 2 standard data structures:
# 1) A SummarizedExperiment will hold primary data from an RNA-Seq pipeline
# e.g. counts, FPKM, TPM, zFPKM, as well as row and column annotation.
#
# 2) Limma output from model fitting will be contained within a SubsettableListOfArrays
# object.  Unlike SummarizedExperiment, this object can hold Fit objects as well as
# multiple types of row and column annotation.
#
# Both SummarizedExperiment and SubsettableListOfArrays object support subsetting the
# all the contained objects appropriately when the main object is subsetted.
# For example:  MySummarizedExperiment[1:100, 1:10] would automatically subset
# the count, FPKM etc matices to the first 100 genes and 10 samples.  The associated
# row and col data is also subsetted appropriately. SubsettableListOfArrays behaves
# similar in this regard.
#
#
# Dealing with multiple contrasts (List of Contrasts dataframes)
#
# Finally, after fitting a model, typically you'll be running multiple contrasts e.g. by
# running make.contrasts and topTable.  It is recommended that your store multiple
# contrasts (topTable output) in a simple named list.  This will facilitate operations using lapply
# to e.g. extract all LogFC columns for a heatmap or filter each contrast in a systematic
# way.
# For example:
#   MyContrastList = list(contrast1 = topTable (....),
#                         contrast2 = topTable (...),
#                         contrast3 = topTable (...)
#                         )
# Then something like:
#   MyLogFCMatrix = MyContrastList %>% lapply(`[[`, "LogFC") %>% do.call(what=cbind)
#
# will aggregate all the logFC columns into one matrix (you'll probably want to redefine
# the column names at this point though)

### DATA STRUCTURES ###

#Need separate SummarizedExperiments for Gene and Transcript level because they have
#different row counts The GeneData list of lists defines the datafiles to load into
#a SummarizedExperiment, the standard name they will be given in the SummarizedExperiment
#object and the type of data contained in that file (rowRanges, colData, Assay or metadata)
#
#Where:
# rowRanges is genes or transcripts annotation
# colData is samples annotation
# Assay is nrow(rowdata) rows by nrow(colData) columns (e.g. counts, FPKM, zFPKM, TPM, etc.)
# The pre-defined GeneData and TranscriptData lists are configured for the Omicsoft pipeline default values.
# You can change the filenames to support other datasources.  But do not change the name or type values.
#
#' GeneData (list of lists)
#'
#' Defines the data files to import into a SummarizedExperiment.
#' Each file is described by a list with defines the fieldname, the type of data
#' and the text file containing the data.  The predefined GeneData item
#' is preconfigured for datafiles from the Omicsoft pipeline.  Modify GeneData
#' file values to use data from a different pipeline.
#'
"GeneData"
GeneData = list(
  Annotation = list(name = "Annotation",
                    file = "RNA-Seq.Count.Annotation.txt",
                    type = "rowRanges"),
  Design = list(name = "Design",
                file = "RNA-Seq.Design.txt",
                type = "colData"),
  Counts = list(name = "Counts",
                file = "RNA-Seq.Count.Table.txt",
                type = "Assay"),
  Level = list(name = "Level",
               level = "Gene",
               type = "metadata"),
  FPKM = list(name = "FPKM",
              file = "RNA-Seq.FPKM.Table.txt",
              type = "Assay"),
  TPM = list(name = "TPM",
              file = "Genes.TPM.Table.txt",
              type = "Assay"),
  QC.Metrics = list(name = "QC.Metrics",
            file = "RNA-Seq.QCMetrics.Table.txt",
            type = "metadata")
)
# uncomment this block to save changes to this datastructure to be built into the package
# x = getwd()
# setwd ("~/R/lib/pkgsrc/DGE.Tools/")
# save(GeneData, file="./data/GeneData.rda")
# setwd(x)

#' TranscriptData (list of lists)
#'
#' Defines the data files to import into a SummarizedExperiment.
#' Each file is described by a list with defines the fieldname, the type of data
#' and the text file containing the data.  The predefined TranscriptData item
#' is preconfigured for datafiles from the Omicsoft pipeline.  Modify TranscriptData
#' file values to use data from a different pipeline.
#'
"TranscriptData"
TranscriptData = list(
  Annotation = list(name = "Annotation",
                    file = "RNA-Seq.Transcript_Count.Annotation.txt",
                    type = "rowRanges"),
  Design = list(name = "Design",
                file = "RNA-Seq.Design.txt",
                type = "colData"),
  Counts = list(name = "Counts",
                file = "RNA-Seq.Transcript_Count.Table.txt",
                type = "Assay"),
  Level = list(name = "Level",
               level = "Transcript",
               type = "metadata"),
  FPKM = list(name = "FPKM",
              file = "RNA-Seq.Transcript_FPKM.Table.txt",
              type = "Assay"),
  TPM = list(name = "TPM",
              file = "Transcripts.TPM.Table.txt",
              type = "Assay"),
  QC.Metrics = list(name = "QC.Metrics",
            file = "RNA-Seq.QCMetrics.Table.txt",
            type = "metadata")
)
# x = getwd()
# setwd ("~/R/lib/pkgsrc/DGE.Tools/")
# save(TranscriptData, file="./data/TranscriptData.rda")
# setwd(x)


###  SubsettableListOrArrays data class
#
#This block defines a subclass of SubsettableListOfArrays.  SubsettableListOfArrays is a
#class object under development by Ryan Thompson at Scripps Research Institute.
#A SubsettableListOfArrays object allows you to define slots for row(gene/Transcript) data
#(e.g. gene annotatiton), column(Sample) data (e.g. Sample metadata) and gene expression
#data with N rows(genes) x M columns (samples).  SummarizedExperiments only allow one
#dataframe each for row and col data.  SubsettableListOfArrays allows you to define
#an arbitrary number of row and col dataframes.  Both SummarizedExperiments and
#SubsettableListOfArrays share the convenient feature that subsetting the main object
#(genes or samples) will properly subset all of the elements contained within.
#

# SLOA <- defineSLOASubclass("SLOA",
#                            I="RowRanges",
#                            IX=c("rowData", "Fit"),  #objects indexed as rowdata
#                            JX=c("colData",  "DesignMatrix"), #objects indexed as col data
#                            IJ=c("RSE", "DGElist", "Elist", "Counts", "FPKM", "Log2FPKM", "Log2CPM", "Residuals"), #data matrices
#                            required=c("rowData", "colData", "Counts")
# )
# I: defines an array of gene data
# IX: defines a dataframe of gene data
#     rowData is for gene annotation
#     Fit will hold a limma fit object
# IJ: define gene x samples data
#     RSE is a SummarizedExperiment (A SummarizedExperiment can be embedded within a SLOA)
#     DGElist, Elist, Counts, and CPM are all intermediate matices that you might
#     want to capture while running a DGE analysis.
# required:  defines which objects must be included in an SLOA
# J:  not used here, but defines an array of sample data
#Now to load data:
#   MyData = SLOA(rowData=MyGeneAnnotation, colData=MySampleAnnotation, Fit = MyFit, Alignment.QC.summary = MyQCSummary, Counts=MyCounts)
#to extract data, reference just like a list:
#   MyCounts = MyData$Counts

### Utility FUNCTIONS ###

### Function tsmsg ###
#' @export
tsmsg <- function(...) {
  # Works like message() but prepends a timestamp
  message(date(), ": ", ...)
}


### Function df2GR ###
#' @import magrittr IRanges GenomicRanges
#' @export
df2GR <- function(df, seqnames=c("seqnames", "chr", "chromosome"), start="start", end="end", strand="strand",
                  start.offset=1, end.offset=start.offset) {
  #Convert a Annotation DF to a genomic Ranges object
  #Optional parameters for seqnames, start, end, strand anticipate possible you might have for these
  #fields in your annotation file.  Only need to modify these if your annotation uses differenc colnames
  #for these fields
  #
    #These lines return the colnames used in your datafile.
    seqnames.col <- match(seqnames, tolower(colnames(df))) %>% na.omit %>% .[1]
    start.col <- match(start, tolower(colnames(df))) %>% na.omit %>% .[1]
    end.col <- match(end, tolower(colnames(df))) %>% na.omit %>% .[1]
    strand.col <- match(strand, tolower(colnames(df))) %>% na.omit %>% .[1]
    other.cols <- setdiff(seq_along(colnames(df)), c(seqnames.col, start.col, end.col, strand.col))

  	#make sure start and end are numeric; if not, remove commas and convert to numeric
  	if (is.character(df[[start.col]])) {
  	  df[[start.col]] = gsub(",", "", df[[start.col]]) %>% as.numeric
  	}
    if (is.character(df[[end.col]])) {
      df[[end.col]] = gsub(",", "", df[[end.col]]) %>% as.numeric
    }

    MyRanges = IRanges::IRanges(start=df[[start.col]] - start.offset + 1,
                                end=df[[end.col]]) - end.offset + 1

    gr <- GenomicRanges::GRanges(seqnames=df[[seqnames.col]],
                                 ranges=MyRanges,
                                 strand=df[[strand.col]])

#     gr <- GenomicRanges::GRanges(seqnames=df[[seqnames.col]],
#                   ranges=IRanges::IRanges(start=df[[start.col]] - start.offset + 1,
#                                  end=df[[end.col]]) - end.offset + 1,
#                   strand=df[[strand.col]])

    GenomicRanges::mcols(gr) <- df[other.cols]
    names(gr) <- rownames(df)
    return(gr)
}

### Function Txt2DF ###
#' @export
Txt2DF <- function(filename) {
  #configured to read Omicsoft .txt files correctly capturing GeneIDs as rownames
  if (file.exists(filename)) {
    df = read.table (filename, sep="\t", stringsAsFactors = FALSE,
                     header=TRUE, row.names = 1, comment.char="",
                     quote="", na.strings=c("NA", "."))
    return (df)
  } else {
    warning (paste ("Warning: File = ", filename, "not found."))
    return (-1)
  }
}
#' Function  eBayes_autoprop
#'
#' If you use the following function in place of eBayes for single contrasts,
#' you can actually use the B-statistic from the results of topTable.
#' The B-statistic is the log-odds of differential expression, so a zero
#' means a 50-50 chance of DE, positive means more likely DE than not,
#' and negative means more likely non-DE. (Note that this statistic uses
#' the natural logarithm, not base 10 or base 2 logs.) However, this
#' statistic is only accurate if you supply a prior estimate of the
#' fraction of genes that are DE. The default proportion is 1%,
#' which should give conservative values as long as at least 1% of
#' genes are DE.
#'
#' The function below uses propTrueNull to determine this prior, which
#' means that the B-statistics should be reliable as long as propTrueNull
#' returns a reasonable-looking value, which you can verify by looking at
#' the p-value distribution.
#'
#' @author Ryan Thompson, \email{rct@@thompsonclan.org}
#' @keywords eBayes, DGE, limma
#'
#' @param ... See ?eBayes
#' @param prop.method defaults to "lfdr"
#'
#' @return SLOA a SubsettableListOfArrays containing data ready for DGE analysis
#'
#' @examples
#' MyFit = eBayes_autoprop (MyFit)
#'
#' @import limma
#'
#' @export
eBayes_autoprop <- function(..., prop.method="lfdr") {
  eb <- eBayes(...)
  ptn <- propTrueNull(eb$p.value, method=prop.method)
  eBayes(..., proportion=1-ptn)
}


#' Function  getLevel
#'
#' Return the Level (Gene or Transcript or Exon) of an RSE object
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RangedSummarizedExperiment
#'
#' @param data Either an RangedSummarizedExperiment, SubsettableListOfArrays or
#'   ContrastList object. A RangedSummarizedExperiment object
#'
#' @return Level of data (Gene, Transcript or Exon)
#'
#' @examples
#' MyLevel = getLevel(RSE)
#' MyLevel <- getLevel(MySLOA)
#' MyLevel <- getLevel(MyContrastData)
#'
#' @import SummarizedExperiment
#'
#' @export
getLevel <- function(data) {
  if (class(data)[[1]] == "RangedSummarizedExperiment"){
    level <- metadata(data)[["Level"]]
  } else {
    level <- data$Level
  }
  return(level)
}
