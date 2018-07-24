# DGE_Tools.R
# Author: JRT 07Feb2017
#
# A set of functions to facilitate differential gene expression analysis.



#Need separate data objects for Gene and Transcript level because they have
#different row counts
#
#The GeneData list of lists defines: the datafiles to load into
#a DGEobj, the standard name they will be given in the SummarizedExperiment
#object and the type of data contained in that file (rowRanges, colData, Assay or metadata)
#
#Where:
# rowRanges is genes or transcripts annotation
# colData is samples annotation
# Assay is nrow(rowdata) rows by nrow(colData) columns (e.g. counts, FPKM, zFPKM, TPM, etc.)
# The pre-defined GeneData and TranscriptData lists are configured for the Omicsoft pipeline default filenames.
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
# "GeneData"
# .GeneData = list(
#   Annotation = list(name = "Annotation",
#                     file = "RNA-Seq.Count.Annotation.txt",
#                     type = "rowRanges"),
#   Design = list(name = "Design",
#                 file = "RNA-Seq.Design.txt",
#                 type = "colData"),
#   Counts = list(name = "Counts",
#                 file = "RNA-Seq.Count.Table.txt",
#                 type = "Assay"),
#   Level = list(name = "Level",
#                level = "Gene",
#                type = "metadata"),
#   # FPKM = list(name = "FPKM",
#   #             file = "RNA-Seq.FPKM.Table.txt",
#   #             type = "Assay"),
#   # TPM = list(name = "TPM",
#   #             file = "Genes.TPM.Table.txt",
#   #             type = "Assay"),
#   QC.Metrics = list(name = "QC.Metrics",
#             file = "RNA-Seq.QCMetrics.Table.txt",
#             type = "metadata")
# )

# uncomment this block to save changes to this datastructure to be built into the package
# x = getwd()
# setwd ("~/R/lib/pkgsrc/DGE.Tools2/")
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
# .TranscriptData = list(
#   Annotation = list(name = "Annotation",
#                     file = "RNA-Seq.Transcript_Count.Annotation.txt",
#                     type = "rowRanges"),
#   Design = list(name = "Design",
#                 file = "RNA-Seq.Design.txt",
#                 type = "colData"),
#   Counts = list(name = "Counts",
#                 file = "RNA-Seq.Transcript_Count.Table.txt",
#                 type = "Assay"),
#   Level = list(name = "Level",
#                level = "Transcript",
#                type = "metadata"),
#   # FPKM = list(name = "FPKM",
#   #             file = "RNA-Seq.Transcript_FPKM.Table.txt",
#   #             type = "Assay"),
#   # TPM = list(name = "TPM",
#   #             file = "Transcripts.TPM.Table.txt",
#   #             type = "Assay"),
#   QC.Metrics = list(name = "QC.Metrics",
#             file = "RNA-Seq.QCMetrics.Table.txt",
#             type = "metadata")
# )
# x = getwd()
# setwd ("~/R/lib/pkgsrc/DGE.Tools2/")
# save(TranscriptData, file="./data/TranscriptData.rda")
# setwd(x)

### Utility FUNCTIONS ###

### Function tsmsg ###
# a timestamped message
#' @export
tsmsg <- function(...) {
  # Works like message() but prepends a timestamp
  message(date(), ": ", ...)
}

#df2GR is obsolete.  Can now use:  gr <- as(geneInfoDF, "GRanges")
### Function df2GR ###
#' @import magrittr IRanges GenomicRanges
df2GR <- function(df, seqnames=c("seqnames", "chr", "chromosome"),
                  start="start", end="end", strand="strand",
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

    GenomicRanges::mcols(gr) <- df[other.cols]
    names(gr) <- rownames(df)
    return(gr)
}

### Function Txt2DF ###
Txt2DF <- function(filename) {
  #configured to read Omicsoft .txt files correctly capturing GeneIDs as rownames
  if (file.exists(filename)) {
    df = read.table (filename, sep="\t", stringsAsFactors = FALSE,
                     header=TRUE, row.names = 1, comment.char="",
                     quote="", na.strings=c("NA", "."),
                     check.names=TRUE)
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


#' Function  yrange
#'
#' extract the Y upper and lower limits from a ggplot2 v3 plot object.
#'
#' @author John Thompson, \email{rct@@thompsonclan.org}
#' @keywords ggplot, ranges, limits
#'
#' @param g A ggplot plot object (ggplot2 v3 or higher)
#'
#' @return A vector length 2
#'
#' @examples
#' myYrange = yrange (myggplot)
#'
#' @import ggplot2
#'
#' @export
#https://gist.github.com/tomhopper/9076152  ranges for ggplot2 v2
yrange <- function(my.ggp){ #pass a ggplot object, return yrange
  # ggplot2 v2 solution:
  # ggplot_build(my.ggp)$layout$panel_ranges[[1]]$y.range
  # ggplot2 v3 solution:
  ggplot_build(my.ggp)$layout$panel_params[[1]]$y.range
}

#' Function  xrange
#'
#' extract the X upper and lower limits from a ggplot2 v3 plot object.
#'
#' @author John Thompson, \email{rct@@thompsonclan.org}
#' @keywords ggplot, ranges, limits
#'
#' @param g A ggplot plot object (ggplot2 v3 or higher)
#'
#' @return A vector length 2
#'
#' @examples
#' myYrange = yrange (myggplot)
#'
#' @import ggplot2
#'
#' @export
xrange <- function(my.ggp){
  # ggplot2 v2:
  # ggplot_build(my.ggp)$layout$panel_ranges[[1]]$x.range
  # ggplot2 v3 solution:
  ggplot_build(my.ggp)$layout$panel_params[[1]]$x.range
}

#footnote
addFootnote <- function (my.ggp, footnoteText, footnoteSize, footnoteColor, footnoteJust=1, yoffset=0){
  #add a right justified (by default) footnote at the bottom plot.
  #footnoteJust: value = 0.1; <0.5 is left justified; > 0.5 is right justified; 0.5 is centered
  #yoffset is fraction of y delta to add to yr[1]

  yr <- yrange(my.ggp)
  xr <- xrange(my.ggp)
  yoffset <- yoffset * (yr[2]-yr[1])
  xcoord <- ifelse(footnoteJust < 0.50, xr[1], xr[2])
  if (footnoteJust == 0.5) #special case = center
    xcoord <- xr[1] + ((xr[2]-xr[1])/2)
  my.ggp <- my.ggp +
    annotate("text", label = footnote, x = xcoord, y = yr[1]+yoffset,
             size = footnoteSize,
             colour = footnoteColor,
             hjust=footnoteJust,
             vjust=1)
}

