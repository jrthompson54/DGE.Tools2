### Function runContrasts ###
#' Function  runContrasts
#'
#' Take an SLOA object and a named list of contrasts to perform. The SLOA must
#' contain a Fit object and a DesignMatrix. Returns a list of data objects
#' including topTable and topTreat output, ContrastMatrix, Fit.Contrasts, gene
#' and sample annotation (rowData and colData).
#'
#' The contrasts are composed of column names from the DesignMatrix of the SLOA
#' and put into list form to start the analysis. Each contrast is named to give
#' it a short pneumonic name.
#'
#' Example ContrastList \cr
#'
#' ContrastList = list( \cr
#'    T1 = "treatment1 - control", \cr
#'    T2 = "treatment2 - control" \cr
#' ) \cr
#'
#'where treatment1, treatment2 and control are columns in the DesignMatrix
#'
#'The returned ContrastAnalysis list contains the followin objects:
#'\itemize{
#'    \item{"ContrastMatrix"} {a matrix}
#'    \item{"Fit.Contrasts"} {a Fit object}
#'    \item{"TopTableList"} {a List of dataframes}
#'    \item{"TopTreatList"} {a List of dataframes}
#'    \item{"SigCountsTable"} {a dataframe summarizing signature gene counts}
#'    \item{"workflowRecord"} {List of reproducible documentation}
#'    \item{"Level"} {Type of data; one of Gene, Transcript or Exon}
#'}
#'
#' See ?extractCol for a convenient means to extract a named column from all
#' the contrasts in a TopTableList or TopTreatList.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords SLOA; contrasts
#'
#' @param sloa An SLOA object containing a Fit object
#' @param ContrastList A named list of contrasts
#' @param runTopTable runs topTable on the specified contrasts (Default = TRUE)
#' @param runTopTreat runs topTreat on the specified contrasts (Default = FALSE)
#' @param FoldChangeThreshold Only applies to TopTreat (Default = 1.5)
#' @param runEBayes Runs eBayes after lmfit; default = TRUE
#' @param robust eBayes robust option (default = TRUE)
#' @param PvalueThreshold Default = 0.01
#' @param FDRthreshold Default = 0.1
#' @param proportion Proportion of genes expected to be differentially expressed
#'   (used by eBayes) (Default = 0.01)
#' @param SigTable Deprecated and ignored: replaced with summarizeSigCounts function
#' @param SigTablePNG Deprecated and ignored: replaced with summarizeSigCounts function
#' @param qvalues Deprecated and ignored.  Qvalue functionality has been removed from
#'   runContrasts and relegated to a separate function see ?runQvalues.#'
#' @return A list containing contrast dataframes and other associated data.frames.
#'
#' @examples
#' #run defaults
#' MyContsratList  = runContrasts (MySLOA, ConstrastList)
#' MyContsratList  = runContrasts (MySLOA, ConstrastList, runTopTable = TRUE
#'        ConstrastType = "TopTreat", FoldChangeThreshold = 1.25)
#'
#' @import S4Vectors SummarizedExperiment zFPKM dplyr limma gridExtra
#'
#' @export
runContrasts <- function(sloa, ContrastList,
									runTopTable = TRUE,
									runTopTreat = FALSE,
                         	FoldChangeThreshold = 1.5,
                         	runEBayes = TRUE,
                         	robust = TRUE,
                         	qvalues = FALSE,
                         	PvalueThreshold=0.01,
                         	FDRthreshold=0.1,
                         	proportion=0.01,
									SigTable=NULL,  #Deprecated and nonfunctional (kept only for backward compatibility)
									SigTablePNG = NULL) {

  FVersion = "runContrasts : 07Jan2016"

  if (!exists("SLOA")) {
    stop("You need to source SubsettableListOfArrays.R before executing runContrasts")
  }
  #argument checks
  if (runTopTable == FALSE & runTopTreat == FALSE){
    stop ("Error: Both runTopTable and runTopTreast set to False; nothing to do")
  }
  if (FoldChangeThreshold <= 0) {
    stop ("Error: FoldChangeThreshold must be a positive value")
  }
  if (!class(ContrastList) == "list") {
    stop ("ContrastList must be class \"list\".")
  }
  if (is.null(names(ContrastList))){
    print ("ContrastList must be a named list but is missing names.")
    print ("Example:")
    print ("   ContrastList = list(")
    print ("     name1 = \"treatment1 - control\",")
    print ("     name2 = \"treatment2 - control\"")
    print ("   )")
    stop ()
  }
  if (qvalues) {
    warning("warning: qvalues is deprecated. Now use runQvalue on the runContrast output if you want to add qvalues")
  }

  #run the contrast fit
  ContrastMatrix = makeContrasts (contrasts=ContrastList, levels=sloa$DesignMatrix)
  MyFit.Contrasts = contrasts.fit(sloa$Fit, ContrastMatrix)
  #run eBayes
  if (runEBayes) {
    MyFit.Contrasts = eBayes(MyFit.Contrasts, robust=robust, proportion=proportion)
    MyFit.Contrasts.treat = treat(MyFit.Contrasts, lfc=log2(FoldChangeThreshold),
    											robust=robust)
  }

  #run topTable on each contrast and add each DF to a list

  if (runTopTable){
    #Run topTable via lapply to generate a bunch of contrasts.
    MyCoef = 1:length(ContrastList) %>% as.list
    MyTopTableList = lapply (MyCoef, function(x) (topTable(MyFit.Contrasts, coef=x,
                              confint=T, number=Inf, p.value=1, sort.by="none")))
    #transfer the contrast names
    names(MyTopTableList) = names(ContrastList)

    #construct a vector of counts that meet thresholds of significance
    # Pval01count = sapply(MyTopTableList, . %$% P.Value %>% is_weakly_less_than(PvalueThreshold) %>% sum)
    # FDR10count = sapply(MyTopTableList, . %$% adj.P.Val %>% is_weakly_less_than(FDRthreshold) %>% sum)
  }

  if (runTopTreat) {
    #Run topTreat via lapply to generate a bunch of contrasts.
    LFC = log2(FoldChangeThreshold)
    MyCoef = 1:length(ContrastList) %>% as.list
    MyTopTreatList = lapply (MyCoef, function(x) (topTreat(MyFit.Contrasts.treat, coef=x,
                                  confint=T, lfc=LFC, number=Inf, p.value=1, sort.by="none")))
    #transfer the contrast names
    names(MyTopTreatList) = names(ContrastList)

    #construct a vector of counts that meet thresholds of significance
    # Pval01countTreat = sapply(MyTopTreatList, . %$% P.Value %>% is_weakly_less_than(PvalueThreshold) %>% sum)
    # FDR10countTreat = sapply(MyTopTreatList, . %$% adj.P.Val %>% is_weakly_less_than(FDRthreshold) %>% sum)
  }

#   #Optionally Provide a data table of signature Counts
#   #
#   #### Add Qvalues to the summary table ###
#   if (SigTable == TRUE) {
#
#     PFCcol = paste("P<", PvalueThreshold, " & ", FoldChangeThreshold, "FC", sep="")
#     FDRFCcol = paste((FDRthreshold*100), "% FDR & ", FoldChangeThreshold, "FC", sep="")
#     PVcol = paste("P<", PvalueThreshold, sep="")
#     FDRcol = paste ((FDRthreshold*100), "% FDR", sep="")
#     LFDRcol = paste ((FDRthreshold*100), "% LFDR", sep="")
#     #put the counts into a dataframe and print
#     library(gridExtra)
#     if (runTopTable == TRUE & runTopTreat == TRUE & qvalues == TRUE){
#       SigCounts = data.frame(cbind(Pval01count, FDR10count, LFDR10count, Pval01countTreat, FDR10countTreat))
#       colnames(SigCounts) = c(PVcol, FDRcol, LFDRcol, PFCcol, FDRFCcol)
#     } else if (runTopTable == TRUE & runTopTreat == TRUE){
#       SigCounts = data.frame(cbind(Pval01count, FDR10count, Pval01countTreat, FDR10countTreat))
#       colnames(SigCounts) = c(PVcol, FDRcol, PFCcol, FDRFCcol)
#     } else if (runTopTable == TRUE & qvalues == TRUE){
#       SigCounts = data.frame(cbind(Pval01count, FDR10count, LFDR10count))
#       colnames(SigCounts) = c(PVcol, FDRcol, LFDRcol)
#     } else if (runTopTable == TRUE){
#       SigCounts = data.frame(cbind(Pval01count, FDR10count))
#       colnames(SigCounts) = c(PVcol, FDRcol)
#     } else if (runTopTreat == TRUE){
#       SigCounts = data.frame(cbind(Pval01countTreat, FDR10countTreat))
#       colnames(SigCounts) = c(PFCcol, FDRFCcol)
#     }
#
#     plot.new()
#     grid.table(SigCounts)
#
#     if (!is.null(SigTablePNG)) {
#       png(file=SigTablePNG, width=8,height=6, units = 'in', res = 300)
#       plot.new()
#       grid.table(SigCounts)
#       invisible ( dev.off() )
#     }
#   }

  workflowRecord = list()
  workflowRecord$Session.Info = sessionInfo()
  workflowRecord$runContrasts = paste (date(), " : v", packageVersion("DGE.Tools"), sep="")
  workflowRecord$Date = date()
  workflowRecord$RVersion = R.version.string

  #Output
  ContrastAnalysis = list(
    ContrastMatrix = ContrastMatrix, #matrix
    Fit.Contrasts = MyFit.Contrasts, #Fit object
    colData = sloa$colData,  #transfer sample annotation
    rowData = sloa$rowData   #transfer gene annotation
  )
  if (runEBayes == TRUE && runTopTreat == TRUE){
  		ContrastAnalysis$Fit.Contrasts.treat = MyFit.Contrasts.treat
  }

  if (runTopTable) {
    ContrastAnalysis$TopTableList = MyTopTableList   #List of DF
  }
  if (runTopTreat) {
    ContrastAnalysis$TopTreatList = MyTopTreatList   #List of DF
  }
  #ContrastAnalysis$SigCountsSummary = SigCounts    #DF
  ContrastAnalysis$workflowRecord = workflowRecord                #List
  ContrastAnalysis$Level = sloa$Level

  return(ContrastAnalysis)
}


