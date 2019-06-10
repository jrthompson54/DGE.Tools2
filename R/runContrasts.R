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
#'   to "fitName_cf".  Change only if you need to create 2 or more contrast sets from the same fit.
#' @param runTopTable runs topTable on the specified contrasts (Default = TRUE)
#' @param runTopTreat runs topTreat on the specified contrasts (Default = FALSE)
#' @param FoldChangeThreshold Only applies to TopTreat (Default = 1.5)
#' @param runEBayes Runs eBayes after contrast.fit (default = TRUE)
#' @param robust eBayes robust option (default = TRUE)
#' @param PvalueThreshold Default = 0.01
#' @param FDRthreshold Default = 0.1
#' @param proportion Proportion of genes expected to be differentially expressed
#'   (used by eBayes) (Default = 0.01)
#' @param Qvalue Set TRUE to include Qvalues in topTable output (Default=FALSE)
#' @param IHW Set TRUE to add FDR values from the IHW package (Defulat=FALSE)
#' @param verbose Set TRUE to print some information during processing (Default=FALSE)
#' @return The DGEobj with contrast fits and topTable/topTreat dataframes added.
#'
#' @examples
#' #run defaults
#' myDgeObj  = runContrasts (myDgeObj, myFitName, ConstrastList)
#' myDgeObj  = runContrasts (myDgeObj, myFitName, ConstrastList, runTopTable = TRUE
#'        runTopTreat=TRUE, FoldChangeThreshold = 1.25)
#'
#' @import magrittr
#' @importFrom limma contrasts.fit eBayes makeContrasts topTable topTreat treat
#' @importFrom DGEobj addItem getItem
#' @importFrom assertthat assert_that
#' @importFrom stringr str_c
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
                         proportion=0.01,
                         Qvalue=FALSE,
                         IHW=FALSE,
                         verbose=FALSE) {

  assertthat::assert_that (!missing(dgeObj),
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
  designMatrix <- try({getItem(dgeObj, designMatrixName)}, silent=TRUE)
  if (class(designMatrix) == "try-error")
    stop(paste("Couldn't find", designMatrixName, "in dgeObj.", sep=" "))

  fitName <- paste(designMatrixName, "_fit", sep="")
  fit <- try((getItem(dgeObj, fitName)), silent=TRUE)
  if (class(fit) == "try-error") {
    stop(paste(fitName, "not found in dgeObj", sep=" "))
  }

  #run the contrast fit
  ContrastMatrix <- limma::makeContrasts (contrasts=contrastList, levels=designMatrix)
  MyFit.Contrasts <- limma::contrasts.fit(fit, ContrastMatrix)

  #run eBayes
  if (runEBayes) {
    if (verbose == TRUE) tsmsg(stringr::str_c("running EBayes: proportion = ", proportion))
    MyFit.Contrasts = limma::eBayes(MyFit.Contrasts, robust=robust, proportion=proportion)
    MyFit.Contrasts.treat = limma::treat(MyFit.Contrasts, lfc=log2(FoldChangeThreshold),
                                  robust=robust)
  }

  #run topTable on each contrast and add each DF to a list

  if (runTopTable == TRUE){
    if (verbose == TRUE) tsmsg("Running topTable...")
    #Run topTable via lapply to generate a bunch of contrasts.
    MyCoef = 1:length(contrastList) %>% as.list
    TopTableList = lapply (MyCoef, function(x) (limma::topTable(MyFit.Contrasts, coef=x,
                              confint=T, number=Inf, p.value=1, sort.by="none")))

    #transfer the contrast names
    names(TopTableList) = names(contrastList)

    if (Qvalue == TRUE){
        TopTableList <- runQvalue(TopTableList)
    }
    if (IHW == TRUE){
        IHW_result <- runIHW(TopTableList)
        TopTableList <- IHW_result[[1]]
    }
  }

  if (runTopTreat == TRUE) {
    if (verbose == TRUE) tsmsg("Running topTreat...")
    #Run topTreat via lapply to generate a bunch of contrasts.
    LFC = log2(FoldChangeThreshold)
    MyCoef = 1:length(contrastList) %>% as.list
    TopTreatList = lapply (MyCoef, function(x) (limma::topTreat(MyFit.Contrasts.treat, coef=x,
                                  confint=T, lfc=LFC, number=Inf, p.value=1, sort.by="none")))
    #transfer the contrast names
    names(TopTreatList) = names(contrastList)
  }

#capture the contrastMatrix, MyFit.Contrasts and the contrast DFs

 #put results into the dgeobj
 #capture the contrast matrix
  dgeObj <- DGEobj::addItem(dgeObj, item=ContrastMatrix,
                            itemName=paste(contrastSetName, "_cm", sep=""),
                            itemType="contrastMatrix", funArgs=funArgs,
                            parent=fitName)


 if (runTopTable){
    #add the contrast fit
   dgeObj <- DGEobj::addItem(dgeObj, item=MyFit.Contrasts,
                             itemName=paste(contrastSetName, "_cf", sep=""),
                             itemType="contrast_fit",
                             funArgs=funArgs,
                             parent=fitName)

    #add the topTable DFs
    listNames <- names(TopTableList)
    for (i in 1:length(TopTableList))
      dgeObj <- DGEobj::addItem(dgeObj, item=TopTableList[[i]],
                                itemName=listNames[[i]],
                                itemType="topTable", funArgs=funArgs,
                                parent=paste(contrastSetName, "_cf", sep=""))
 }

 if (runTopTreat){
   dgeObj <- DGEobj::addItem(dgeObj, item=MyFit.Contrasts.treat,
                             itemName=paste(contrastSetName, "_cft", sep=""),
                             itemType="contrast_fit_treat",
                             funArgs=funArgs,
                             parent=fitName)

    listNames <- names(TopTreatList)
    for (i in 1:length(TopTreatList))
      dgeObj <- DGEobj::addItem(dgeObj, item=TopTreatList[[i]],
                                itemName=paste(listNames[i], "_treat", sep=""),
                                itemType="topTreat", funArgs=funArgs,
                                parent=paste(contrastSetName, "_cft", sep=""))
 }

 return(dgeObj)
}


