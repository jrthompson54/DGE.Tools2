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
                              customAttr){
    
    #get the data i
    countData <- read.table (counts, sep="\t", stringsAsFactors = FALSE,
                             header=TRUE, row.names = 1, comment.char="",
                             quote="", na.strings=c("NA", "."))
    
    seqData <- read.table (seqAnnotation, sep="\t", stringsAsFactors = FALSE,
                             header=TRUE, row.names = 1, comment.char="",
                             quote="", na.strings=c("NA", ".")) %>% 
                             as.matrix
    
    designData <- read.table (design, sep="\t", stringsAsFactors = FALSE,
                             header=TRUE, row.names = 1, comment.char="",
                             quote="", na.strings=c("NA", "."))
    
    #build the DgeObj
    if (missing(customAttr))
        DgeObj <- initDGEobj(countData, designData, seqData, level)
    else {
        if (class(customAttr)[[1]] == "list")
            DgeObj <- initDGEobj(countData, designData, seqData, level, 
                                 customAttr=customAttr)
        else
            stop ("customAttr must be a named list of attribute/value pairs")
    
    }
    return(DgeObj)
}

  
