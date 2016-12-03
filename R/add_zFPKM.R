### Function Add_zFPKM ###
#' Function  Add_zFPKM
#'
#' Add a zFPKM assay to the RSE.  Optionally, plot the zFPKM distributions and fits.
#'
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RangedSummarizedExperiment, zFPKM
#'
#' @param RSE A RangedSummarizedExperiment (required)
#' @param outputPath location for the plot file (default = currentDir)
#' @param plotFile Name for image of sample zFPKM distribution fits (default =
#' "zFPKM.png"; set to NULL to disable plot)
#' @param facetTitles default = TRUE. Set to FALSE if facet titles too cluttered.
#'
#' @return A RangedSummarizedExperiment with zFPKM assay added
#'
#' @examples
#' MyRSE = add_zFPKM (MyRSE)
#' MyRSE = add_zFPKM (MyRSE, plotFile = "Gene.zFPKM.PNG")
#' MyRSE = add_zFPKM (MyRSE, plotFile = "Gene.zFPKM.PNG", outputPath="./plots")
#' MyRSE = add_zFPKM (MyRSE, plotFile = NULL)
#' 
#' @import S4Vectors SummarizedExperiment zFPKM magrittr
#'
#' @export
add_zFPKM <- function (RSE, OutputPath=".", 
                       PlotFile = "zFPKM.PNG", 
                       facetTitles=FALSE){
  level <- (extract from metadata slot)
  Warn if level <> gene 
  
  MyRowData <- 
      
  #Calculate RawLog2TPM, and RawLog2FPKM and add to MyAssays if ExonLength present
  if (!is.null(MyRowData$ExonLength)) {


    #we'll use the imported fpkm table to calculate zFPKM
      MyAssays$zFPKM = zFPKMTransformDF (data.frame(MyAssays$FPKM), PlotDir = OutputPath,
                                         PlotFile=plotFile, 
                                         FacetTitles=facetTitles)
      tsmsg("zFPKM calculated from FPKM data and added to assays.")
  }
***Get workflowRecord from RSE

  #Capture reproducible info in metadata
  workflowRecord$Session.Info = sessionInfo()
  workflowRecord$add_zFPKM = paste (date(), " : v", packageVersion("DGE.Tools2"), sep="")
  workflowRecord$RVersion = R.version.string
  workflowRecord$Date = date()

  Mymetadata$workflowRecord = workflowRecord
  
  tsmsg("Appending zFPKM RSE Object")

  ***put zFPKM into assays and metadata back into metadata slot

    return(RSE)

}
