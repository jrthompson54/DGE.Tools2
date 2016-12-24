### Function Build_RSE ###
#' Function  Build_RSE
#'
#' Build a RangedSummarizedExperiment (RSE) from count data.  The
#' RSE is returned.
#'
#' User provides separate tab-delimited text files for each row annotation,
#' column anotation and "assays" (matrices of RxC e.g. counts, FPKM, etc).
#'
#' GeneData is a predefined data structure (list of lists) that specifies the
#' datafiles to import.  GeneData is pre-configured for the standarized Omicsoft
#' project output.  Modify as necessary to accommodate data from other
#' pipelines.  There is also a "TranscriptData" pre-defined definition.  Passing
#' TranscriptData as the GeneData parameter will produce a
#' Transcript-level RSE.  Note a single RSE cannot contain both Gene and Transcript
#' data because of the mismatch in gene/transcript rowcounts.
#'
#' Data Requirements:
#'
#' Minimally, Sample (column) annotation, gene or transcript (row) annotation and
#' counts are required to build an RSE object.
#'
#' Column 1 in annotation and assay data files should be named GeneID or
#' TranscriptID and the ID order in all files should be the same.
#'
#' Similarly, the row order of the Column annotation should match the
#' column order in the assay data files.
#'
#' The Gene/Transcript annotation must include chromosome position data (chr,
#' start, end, strand).  FPKM and TPM will be imported if available or
#' calculated and added to the RSE if an
#' ExonLength field is provided in the row annotation.  Alignment QC data
#' will also be imported, if present, and placed in the metadata slot.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID, SummarizedExperiment, Omicsoft
#'
#' @param GeneData A predefined List of Lists defining data names, file names, and file types.
#' @param OutputPath  Path 
#'
#' @return A RangedSummarizedExperiment
#'
#' @examples
#' 
#'     #For Omicsoft data in same folder
#'     MyRSE <- Build_RSE ()
#'     
#'     #User supplied files and modified .Genedata list accordingly to point
#'     #to their files
#'     MyRSE <- Build_RSE (MyGeneData, OutputPath = "./rdata")
#'
#' @import S4Vectors SummarizedExperiment dplyr magrittr
#'
#' @export
Build_RSE <- function (.GeneData, OutputPath="."){

  MyAssays = list()
  Mymetadata = list()

  #check for level definition
  if (is.null(GeneData$Level)){
    Mymetadata$Level = "Unk"
    stop ("GeneData$Level$level must have be set to Gene or Transcript")
  }
  tsmsg("Importing text files ...")
  #read tabbed text files specified in GeneData into dataframes
  for (i in 1:length(GeneData)) {

    #check for Level flag first,  then read files
    if (GeneData[[i]]$name == "Level") {
      Mymetadata[[GeneData[[i]]$name]] <- GeneData[[i]]$level

    } else if (file.exists(GeneData[[i]]$file)) {
      if  (GeneData[[i]]$type == "Assay") {
        MyAssays[[GeneData[[i]]$name]] <- as.matrix (Txt2DF (GeneData[[i]]$file))
        #10Aug2016 need to use make.names to avoid invalid rownames
        colnames(MyAssays[[GeneData[[i]]$name]]) %<>% make.names
        tsmsg(paste(GeneData[[i]]$name, "added to assays."))
      }else if  (GeneData[[i]]$type == "rowRanges") {
        MyRowData = Txt2DF (GeneData[[i]]$file)
        tsmsg(paste(GeneData[[i]]$name, "added to rowRanges."))
      }else if  (GeneData[[i]]$type == "colData") {
        MyColData = Txt2DF (GeneData[[i]]$file)
        #10Aug2016 need to use make.names to match assay rownames
        rownames(MyColData) %<>% make.names
        tsmsg(paste(GeneData[[i]]$name, "added to colData."))
      }else if  (GeneData[[i]]$type == "metadata") {
        #metadataCount = metadataCount + 1
        Mymetadata[[GeneData[[i]]$name]] = Txt2DF (GeneData[[i]]$file)
        tsmsg(paste(GeneData[[i]]$name, "added to metadata."))
      }
    } else if (GeneData[[i]]$type %in% c("rowRanges", "colData")) {
      stop (sprintf ("%s is required.  Please add the \"%s\" file and rerun.",
                     GeneData[[i]]$type, GeneData[[i]]$file))
    } else {
      warning(paste("Warning: File", GeneData[[i]]$file, "Not found! Skipped."))
    }
  }

  # Reality Checks

  # Gene/transcript Counts Row order should match row Annotation Gene/TranscriptID.
  # compare rownames(Counts) to MyRowData$GeneID
  if (!is.null(MyAssays$Counts)) {
    #check id order Gene Level
    if (tolower(GeneData$Level$level) == "gene") {
      if (sum(!rownames(MyAssays$Counts) == MyRowData$GeneID) > 0){
        #at least one GeneID did not match
        stop ("GeneIDs in Counts do not match GeneIDs in Row Annotation (check sort order maybe?)")
      }
    }
    #browser()
    #check id match for Transcript level
    if (tolower(GeneData$Level$level) == "transcript") {
      if (sum(!rownames(MyAssays$Counts) == rownames(MyRowData)) > 0){
        #at least one GeneID did not match
        stop ("TranscriptIDs in Counts do not match TranscriptIDs in Row Annotation (check sort order maybe?)")
      }
    }

  } else { #no counts found
    stop ("Count data is required but missing.")
  }

  if (!exists("MyRowData")) {
    stop("Error: Row metadata (gene annotation) required but not found!")
  }
  if (!exists("MyColData")) {
    stop("Error: Column metadata (Sample annotation) required but not found!")
  }

  #Calculate some summary stats for counts if Counts is present and add columns to MyRowData
  if (!is.null(MyAssays$Counts)) { #test if $Counts exists

    CountStats = data.frame(row.names= row.names(MyAssays$Counts))
    CountStats$ID = row.names(MyAssays$Counts)
    CountStats$countMin = apply(MyAssays$Counts, 1, min)
    CountStats$countMax = apply(MyAssays$Counts, 1, max)
    CountStats$countMean = apply(MyAssays$Counts, 1, mean)
    CountStats$countMedian = apply(MyAssays$Counts, 1, median)
    CountStats$countVariance = apply(MyAssays$Counts, 1, var)

    #merge CountStats into rowRanges dataframe
    MyRowData$ID = rownames(MyRowData)
    MyRowData = dplyr::left_join(MyRowData, CountStats)
    rownames(MyRowData) = MyRowData$ID
    #MyRowData$GeneID = NULL
    tsmsg("Added CountStats (min, max, mean counts) to rowRanges.")
  }

  #Calculate QC.Metrics.Summary if QC.Metrics is present and add to Mymetadata
  if (!is.null(Mymetadata$QC.Metrics)) {
    QC.Metrics.Summary = data.frame(row.names= row.names(Mymetadata$QC.Metrics))
    QC.Metrics.Summary$min = apply(Mymetadata$QC.Metrics, 1, min)
    QC.Metrics.Summary$max = apply(Mymetadata$QC.Metrics, 1, max)
    QC.Metrics.Summary$mean = apply(Mymetadata$QC.Metrics, 1, mean)
    QC.Metrics.Summary$mean = apply(Mymetadata$QC.Metrics, 1, median)
    Mymetadata$QC.Metrics.Summary = QC.Metrics.Summary
    tsmsg("Added QC.Metrics.Summary to metadata.")
  }
  #A QC.Metrics summary was added to metadata in the for loop above. The QC.Metrics
  #metadata dataframe has rows = metrics and cols = samples (as it comes out of
  #Omicsoft). 
  
  #Capture reproducible info in metadata
  workflowRecord = list()
  workflowRecord$Session.Info = sessionInfo()
  workflowRecord$RVersion = R.version.string
  workflowRecord$DGE.Tools.Version = packageVersion("DGE.Tools2")
  workflowRecord$Date = date()

  Mymetadata$workflowRecord = workflowRecord

  tsmsg("Building RSE Object")


  RSE = SummarizedExperiment(assays = MyAssays,
                               rowRanges = df2GR(MyRowData),
                               colData = S4Vectors::DataFrame(MyColData),
                               metadata = S4Vectors::SimpleList(Mymetadata))

  return(RSE)

}
