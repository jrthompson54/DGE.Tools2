### Function runContrasts ###
#' Function  runContrasts
#'
#' Take a DGEObj and a named list of contrasts to build. The DGEobj must
#' contain a limma Fit object and associated DesignMatrix. Returns the DGEobj with
#' contrast fit(s), contrast matrix and topTable/topTreat dataframes added.  
#'
#' The contrastList is a named list composed of column names from the DesignMatrix 
#' of the DGEobj.  Each contrast is named to give
#' it a short recognizable name.
#'
#' Example contrastList \cr
#'
#' contrastList = list( \cr
#'    T1 = "treatment1 - control", \cr
#'    T2 = "treatment2 - control" \cr
#' ) \cr
#'
#'where treatment1, treatment2 and control are columns in the DesignMatrix
#'
#'The returned ContrastAnalysis list contains the following objects:
#'\itemize{
#'    \item{"ContrastMatrix"} {a matrix}
#'    \item{"Fit.Contrasts"} {a Fit object}
#'    \item{"TopTableList"} {a List of dataframes}
#'    \item{"TopTreatList"} {a List of dataframes}
#'}
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq; DGE; contrasts; DGEobj
#'
#' @param dgeObj A DGEobj object containing a Fit object and design matrix (required)
#' @param designMatrixName The name of the design matrix within dgeObj to use for 
#'    contrast analysis  (required)
#' @param contrastList A named list of contrasts  (required)
#' @param contrastSetName Name for the set of contrasts specified in contrastList.  Defaults
#'   to "fitName_cf".  Change if you need to create 2 or more contrast sets from the same fit.
#' @param runTopTable runs topTable on the specified contrasts (Default = TRUE)
#' @param runTopTreat runs topTreat on the specified contrasts (Default = FALSE)
#' @param FoldChangeThreshold Only applies to TopTreat (Default = 1.5)
#' @param runEBayes Runs eBayes after lmfit (default = TRUE)
#' @param robust eBayes robust option (default = TRUE)
#' @param PvalueThreshold Default = 0.01
#' @param FDRthreshold Default = 0.1
#' @param proportion Proportion of genes expected to be differentially expressed
#'   (used by eBayes) (Default = 0.01)
#' @return The DGEobj with contrast fits and topTable/topTreat dataframes added.
#'
#' @examples
#' #run defaults
#' myDgeObj  = runContrasts (myDgeObj, myFitName, ConstrastList)
#' myDgeObj  = runContrasts (myDgeObj, myFitName, ConstrastList, runTopTable = TRUE
#'        runTopTreat=TRUE, FoldChangeThreshold = 1.25)
#'
#' @import DGEobj dplyr limma gridExtra magrittr
#'
#' @export
runContrasts <- function(dgeObj, designMatrixName, 
                         contrastList,
                         contrastSetName = fitName,
						runTopTable = TRUE,
						runTopTreat = FALSE,
                        FoldChangeThreshold = 1.5,
						PvalueThreshold=0.01,
						FDRthreshold=0.1,
                        runEBayes = TRUE,
                        robust = TRUE,
                        proportion=0.01) {

  assert_that (!missing(dgeObj),
               !missing(designMatrixName),
               !missing(contrastList),
               class(dgeObj)[[1]] == "DGEobj",
               class(contrastList)[[1]] == "list",
               FoldChangeThreshold >=0,
               !is.null(names(contrastList)),
               PvalueThreshold > 0 & PvalueThreshold <= 1,
               !(runTopTable == FALSE & runTopTreat == FALSE)
               )
    
  funArgs <- match.call() #capture arguments
 
  #need to retrieve designMatrix
  result <- try({designMatrix <- getItem(dgeObj, designMatrixName)}, silent=TRUE)
  if (class(result) == "try-error")
      stop(paste("Couldn't find", designMatrixName, "in dgeObj.", sep=" "))
  
  fitName <- paste(designMatrixName, "_fit", sep="")
  result <- try((fit <- getItem(dgeObj, fitName)), silent=TRUE)
  if (class(result) == "try-error") {
      stop(paste(fitName, "not fount in dgeObj", sep=" "))
  }

  #run the contrast fit
  ContrastMatrix <- makeContrasts (contrasts=contrastList, levels=designMatrix)
  MyFit.Contrasts <- contrasts.fit(fit, ContrastMatrix)

  #run eBayes
  if (runEBayes) {
    MyFit.Contrasts = eBayes(MyFit.Contrasts, robust=robust, proportion=proportion)
    MyFit.Contrasts.treat = treat(MyFit.Contrasts, lfc=log2(FoldChangeThreshold),
    											robust=robust)
  }
  
  #run topTable on each contrast and add each DF to a list

  if (runTopTable == TRUE){
    #Run topTable via lapply to generate a bunch of contrasts.
    MyCoef = 1:length(contrastList) %>% as.list
    TopTableList = lapply (MyCoef, function(x) (topTable(MyFit.Contrasts, coef=x,
                              confint=T, number=Inf, p.value=1, sort.by="none")))
    #transfer the contrast names
    names(TopTableList) = names(contrastList)
  }

  if (runTopTreat == TRUE) {
    #Run topTreat via lapply to generate a bunch of contrasts.
    LFC = log2(FoldChangeThreshold)
    MyCoef = 1:length(contrastList) %>% as.list
    TopTreatList = lapply (MyCoef, function(x) (topTreat(MyFit.Contrasts.treat, coef=x,
                                  confint=T, lfc=LFC, number=Inf, p.value=1, sort.by="none")))
    #transfer the contrast names
    names(TopTreatList) = names(contrastList)
  }

#capture the contrastMatrix, MyFit.Contrasts and the contrast DFs

 #capture the contrast matrix
 dgeObj <- addItem(dgeObj, item=ContrastMatrix, 
                   itemName=paste(contrastSetName, "_cm", sep=""),
                   itemType="contrastMatrix", funArgs=funArgs,
                   parent=fitName)
 
 if (runTopTable){
    #add the contrast fit
    dgeObj <- addItem(dgeObj, item=MyFit.Contrasts,
                      itemName=paste(contrastSetName, "_cf", sep=""),
                      itemType="contrast_fit",
                      funArgs=funArgs,
                      parent=fitName)
    
    #add the topTable DFs
    listNames <- names(TopTableList)
    for (i in 1:length(TopTableList))
        dgeObj <- addItem(dgeObj, item=TopTableList[[i]], 
                          itemName=listNames[[i]], 
                          itemType="topTable", funArgs=funArgs,
                          parent=paste(contrastSetName, "_cf", sep=""))
 }
  
 if (runTopTreat){
    dgeObj <- addItem(dgeObj, item=MyFit.Contrasts.treat,
                       itemName=paste(contrastSetName, "_cft", sep=""),
                       itemType="contrast_fit_treat",
                       funArgs=funArgs,
                       parent=fitName) 
    
    listNames <- names(TopTreatList)
    for (i in 1:length(TopTreatList))
        dgeObj <- addItem(dgeObj, item=TopTreatList[[i]], 
                          itemName=paste(listNames[i], "_treat", sep=""), 
                          itemType="topTreat", funArgs=funArgs,
                          parent=paste(contrastSetName, "_cft", sep=""))
 }
     
 return(dgeObj)
}


