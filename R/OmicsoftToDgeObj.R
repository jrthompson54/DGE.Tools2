### Function OmicsoftToDgeObj ###
#' Function  OmicsoftToDgeObj
#'
#' Build a DGEobj from Omicsoft output or tabbed text files.
#'
#' User provides separate tab-delimited text files for row (gene) annotation,
#' column (sample) anotation and "assays" (matrices of RxC e.g. counts, FPKM, etc).
#'
#' Two global variables hold the names of the Omicsoft data files; .geneData
#' and .isoformData.
#'
#' Data Requirements:
#'
#' Assay data should have rownames (sequence ids) and column names (sample IDs).
#' Sequence annotation should also have the same rownames as assays.
#'
#' If possible, the sequence annotation should include chromosome position data (chr,
#' start, end, strand).
#'
#' Sample annotation is one row for each column in the count table.
#' rownames(sampleannotation) == colnames(counts).
#'
#' Function DGEobj::annotateDGEobj provides an easier way than the customAttr argument
#' here.  annotateDGEobj reads key=value pairs from a text file to define
#' attributes.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords Omicsoft, DGEObj, RNA-Seq, Data loading
#'
#' @param path File path for the three data files (Default = "./")
#' @param counts A text file name for count data (gene rows by sample columns) (required)
#'  [Default = "RNA-Seq.Count.Table.txt"]
#' @param seqAnnotation  Filename for Gene, isoform or exon level (row) annotation (required)
#'  [Default = "RNA-Seq.Count.Annotation.txt"]
#' @param design Filename for sample annotation with expt factors and other sample-associated
#'     data (required) [Default = RNA-Seq.Design.txt"]
#' @param level One of "gene", "isoform", "exon" (required) [Default = "gene"]
#' @param source Default = "Omicsoft.  Change if your data if from somewhere else.
#' @param customAttr A named list of custom attributes to attach to the DGEobj;
#'    Optional but highly encouraged;  suggestions: list(PID = "20170101-0001",
#'    XpressID = "123456", Genome = "Mouse.B38", GeneMobel = "Ensembl.R84")
#'
#' @return A DGEobj
#'
#' @examples
#'
#' #Defaults set for an omicsoft dataset:
#' MyDgeObj <- OmicsoftToDgeObj (customAttr=list(
#'                                    PID = "20170101-0001",
#'                                    XpressID = "123456",
#'                                    Genome = "Mouse.B38",
#'                                    GeneMobel = "Ensembl.R84")
#'
#' #I have data files from somewhere else
#' MyDgeObj = OmicsoftToDgeObj (MyCounts, MyGeneAnnotation, MyDesign,
#'                                    level = "gene",
#'                                    customAttr = list(
#'                                       PID = "20170101-0001",
#'                                       XpressID = "123456",
#'                                       Genome = "Mouse.B38",
#'                                       GeneMobel = "Ensembl.R84")
#'                              )
#'
#' @importFrom stringr str_c
#' @importFrom DGEobj initDGEobj
#'
#' @export
OmicsoftToDgeObj <- function (counts = "RNA-Seq.Count.Table.txt",
                              seqAnnotation = "RNA-Seq.Count.Annotation.txt",
                              design = "RNA-Seq.Design.txt",
                              level = "gene",
                              source = "Omicsoft",
                              path = ".",
                              customAttr,
                              gz=FALSE){

    message("OmicsoftToDgeObj is deprecated as of 0.9.56. Use textToDgeObj instead.")
    #add support for gzipped files
    if (gz==TRUE){
      gz <- ".gz"
    } else {
      gz <- ""
    }
    counts <- stringr::str_c(counts, gz)
    seqAnnotation <- stringr::str_c(seqAnnotation, gz)
    design <- stringr::str_c(design, gz)


   #change default filenames if not given and level = isoform
    if (tolower(level) == "isoform") {
        if (missing(counts))
            counts <- stringr::str_c("RNA-Seq.Transcript_Count.Table.txt", gz)
        if (missing(seqAnnotation))
            seqAnnotation <- stringr::str_c("RNA-Seq.Transcript_Count.Annotation.txt", gz)
        if (missing(design))
            design <- stringr::str_c("RNA-Seq.Design.txt", gz)
    }

    #get the data
    countData <- Txt2DF(file.path(path, counts))
    seqData <- Txt2DF(file.path(path, seqAnnotation))
    designData <- Txt2DF(file.path(path, design))
    rownames(designData) <- make.names(rownames(designData))


    #add source to customAttr
    if (missing(customAttr)) {
        customAttr <- list(source="Omicsoft")
    } else {
        assertthat::assert_that(class(customAttr)[[1]] == "list")
        customAttr$source <- source
    }

    #Add DGE.Tools2Version info
    customAttr$DGE.Tools2 <- packageVersion("DGE.Tools2")
    customAttr$DGEobj <- packageVersion("DGEobj")

    #build the DgeObj
    DgeObj <- DGEobj::initDGEobj(counts=countData, rowData=seqData,
                        colData=designData, level,
                        customAttr=customAttr)

    return(DgeObj)
}


