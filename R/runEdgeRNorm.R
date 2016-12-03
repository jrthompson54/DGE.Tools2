### Function runEdgeRNorm ###
#' Function  runEdgeRNorm (edgeR normlization )
#'
#' Returns a SubsettableListOfArrays (SLOA) containing the Fit and associated
#' output from the standard RNA-Seq using edgeR TMM normalization
#' followed by limma voom and lmfit.  The SLOA contains the Fit object
#' and various associated data  use "print(SLOA)" to examine the contents.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID
#'
#' @param RSE An RSE object e.g as built by Build_RSE function
#' @param Formula a text string specifying the model formula. The formula must
#'    use column names from the sample annotation in the RSE.DGE
#' @param NormMethod One of "TMM", "RLE", "upperquartile" or "none". Default = "TMM"
#' @param Design Used to override the design (colData) embedded in the RSE. Defaults to NULL
#' @param plotFile Enable a Barplot of the norm.factors produced (Default = "Norm.Factors.PNG")
#'   Set to NULL to disable the plot.
#' @param plotLabels Sample labels for the plot. Length must equal the number of
#'   samples. (Default = NULL;  sample number will be displayed)
#' @param estimateDispersion If true, run edgeR estimateDisp function and adds
#'   dispersion values to the SLOA output.  (Default = True)  This takes a
#'   minute or two ro run and you may want to turn it off when you re-running a
#'   script during development or maybe better, use the cache=TRUE chunk opion
#'   in knitr so save time rerunning blocks that haven't changed.
#'
#' @return SLOA a SubsettableListOfArrays containing data ready for DGE analysis
#'
#' @examples
#' MyFormula <- "~ 0 + DiseaseTreatment"
#' MySLOA <- runEdgeRNorm (RSE, MyFormula)
#'
#' @import edgeR magrittr limma SummarizedExperiment
#'
#' @export
runEdgeRNorm <- function(RSE, Formula, NormMethod="TMM", Design=NULL,
                         plotFile="Norm.Factors.PNG",
                         plotLabels = NULL,
                         estimateDispersion = TRUE){

  #Version Info
  FVersion = "runEdgeRNorm : 29Dec2015"

  if (!exists("SLOA")) {
    stop("You need to source SubsettableListOfArrays.R before executing runEdgeRNorm")
  }

  #Build design matrix from formula and design data.frame
  #Get design from RSE unless explicitly supplied
  if (is.null(Design)) {
    Design = colData(RSE) %>% as.data.frame
  }
  MyDesignMatrix <- model.matrix(as.formula(Formula), Design)
  #remove colons for interaction terms from colnames (issue #9)
  colnames(MyDesignMatrix) <- make.names(colnames(MyDesignMatrix))

  #need to determine if RSE contains Gene.Counts or Transcript.Counts
  #Get counts from Gene.Counts or Transcript.Counts
  AssayList = assays(RSE) %>% as.list  #get a list of assay matrices
  if (!is.null(AssayList$Counts)){
    MyCounts = assay(RSE, "Counts")
    Level = getLevel(RSE)
  } else {
    Level = NULL
    return -1
    stop("Error: No Counts were found in RSE! Can't proceed")
  }

  #convert counts to a matrix
  CountsMatrix = as.matrix(MyCounts)

  #RNA-Seq Pipeline Block ###############################
  #Now ready to run the DGE pipeline
  MyDGElist = CountsMatrix %>%
    DGEList %>%  #edgeR
    calcNormFactors (method = NormMethod)   #edgeR (does what the function says)
  #MyDGElist is a DGElist object containing:
  #[1] "counts"  "samples"

  if (estimateDispersion == TRUE) {
    MyDGElist %<>% estimateDisp(design = MyDesignMatrix, robust = TRUE)
  }
  #added dispersion parameters to the DGElist:
  #tagwise.dispersion, common.dispersion, trended.dispersion, trend.method,
  #AveLogCPM, span, prior.n, prior.df

  workflowRecord = list()
  workflowRecord$Session.Info = sessionInfo()
  workflowRecord$runEdgeRNorm = paste (date(), " : v", packageVersion("DGE.Tools"), sep="")
  workflowRecord$NormMethod = NormMethod
  workflowRecord$Date = date()
  workflowRecord$RVersion = R.version.string

  #RowData = rowRanges(RSE) %>% elementMetadata %>% as.data.frame
  #Modify to work with rowData and rowRanges
  RowData <- rowData(RSE) %>% as.data.frame
  rownames(RowData) <- RowData$ID

  #fix colnames names in the DesignMatrix (convert colons in interaction terms to periods)
  colnames(MyDesignMatrix) <- make.names(colnames(MyDesignMatrix))

  MySLOA = SLOA(#RowRanges = rowRanges(RSE),
                rowData=RowData,
                colData = colData(RSE),
                Counts = assay(RSE, "Counts"),
                DesignMatrix = MyDesignMatrix,
                DGElist = MyDGElist,
                #other (non-indexed items)
                Formula = Formula,
                QC.Metrics.Summary = metadata(RSE)[["QC.Metrics.Summary"]],
                workflowRecord = workflowRecord,
                Level = Level
  )

  if (!is.null(metadata(RSE)[["QC.Metrics.Summary"]])) { #carry this item over if present
  	MySLOA$QC.Metrics.Summary <- metadata(RSE)[["QC.Metrics.Summary"]]
  }

  #plot the Norm factors
  if (!is.null(plotLabels)  && length(plotLabels == ncol(MySLOA))){
    x = plotLabels
    angle = 45
  } else {
    x = 1:ncol(MySLOA)
    angle = 0
  }

  if (!is.null(plotFile)){ #barplot of norm.factors
    df <- data.frame(x = factor(x),
                     Norm.Factors = MyDGElist$samples$norm.factors)
    nfplot <- ggplot (df, aes(x=x, y=Norm.Factors)) +
                geom_bar (stat="identity", color = "dodgerblue4",
                          fill = "dodgerblue3", width = 0.7) +
                geom_hline(yintercept = 1.0, color = "red") +
                xlab("Samples") +
                ylab("Norm Factors") +
                ggtitle ("Normalization Factors") +
                bwTheme(12) +
                theme(axis.text.x = element_text(angle=angle, hjust = 1.0))

    print(nfplot)
    printAndSave (nfplot, filename = plotFile, printPlot=FALSE)

  }

  return(MySLOA)
}
