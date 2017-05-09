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
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords Omicsoft, DGEObj, RNA-Seq, Data loading
#'
#' @param path File path for the three data files (Default = "./")
#' @param counts A matrix or dataframe of R gene by C samples (required)
#'  [Default = "RNA-Seq.Count.Table.txt"]
#' @param seqAnnotation  Gene, isoform or exon level (row) annotation (required)
#'  [Default = "RNA-Seq.Count.Annotation.txt"]
#' @param design Sample annotation with expt factors and other sample-associated
#'     data (required) [Default = RNA-Seq.Design.txt"]
#' @param level One of "gene", "isoform", "exon" (required) [Default = "gene"]
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
#' @import DGEobj magrittr 
#'
#' @export
OmicsoftToDgeObj <- function (counts = "RNA-Seq.Count.Table.txt", 
                              seqAnnotation = "RNA-Seq.Count.Annotation.txt",
                              design = "RNA-Seq.Design.txt",
                              level = "gene",
                              path = ".",
                              customAttr){
    
    #change default filenames if not given and level = isoform
    if (tolower(level) == "isoform") {
        if (missing(counts))
            counts <- "RNA-Seq.Transcript_Count.Table.txt"
        if (missing(seqAnnotation)) 
            seqAnnotation <- "RNA-Seq.Transcript_Count.Annotation.txt"
        if (missing(design))
            design <- "RNA-Seq.Design.txt"
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
        assert_that(class(customAttr)[[1]] == "list")
        customAttr$source <- "Omicsoft"
   }

    #build the DgeObj
    DgeObj <- initDGEobj(counts=countData, rowData=seqData, 
                        colData=designData, level, 
                        customAttr=customAttr)

    return(DgeObj)
}

 
