### Function runIHW ###
#' Function  runIHW
#'
#' This is a wrapper around the independent hypothesis weighting package that
#' takes a contrast list from runContrasts and applies Independent Hypothesis
#' Weighting (IHW) to the TopTableList within the contrast list.  IHW is a
#' method developed by N. Ignatiadis (http://dx.doi.org/10.1101/034330) to weight
#' FDR values based on a covariate (AveExpr in this case).
#'
#' The IHW FDR values are added as additional columns to the topTable dataframes.
#' It is designed to work on dataframes created by topTable and containe in the
#' "contrast list" data structure produced by runContrasts.
#'
#' Note that it is impractical to run IHW on a list of genes less than ~5000 because
#' it breaks the data into bins of 1500 genes for the analysis if bins = 1, IHW converges
#' on the BH FDR value and there's no point.  So we run IHW on the complete set
#' of detected genes from topTable results (not topTreat).
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2, png, bmp, tiff, jpeg, pdf
#'
#' @param contrastList A list object produced by runContrasts.  The list must
#'   contain, either or both of TopTableList and TopTreatList items as provided
#'   by the runContrasts function.
#' @param alpha The alpha value is the FDR level you wish to interogate (range 0-1; default = 0.1)
#' @param FDRthreshold Value used to estimate proportion of differential genes for each contrast.
#' @param ... other arguments are passed directly to the ihw function (see ?ihw).
#'
#' @return The contrastList is returned.  The TopTableList and TopTreatList items will
#' contain additional columns added by the ihw analysis and prefixed with "ihw."
#' The IHW dataframe of results will also be appended to the contrast list.
#'
#' @examples
#' MyContrasts = runIHW(MyContrasts)
#'
#' @import IHW
#'
#' @export
 runIHW <- function(contrastList,
                    alpha=0.1,
                    FDRthreshold = 0.1,
                    ...){

  FVersion <- "07Jan2016"

  #Function def
    getProportion <- function(ttdf, threshold) {
      #get the proportion for one df
      bhfdrproportion <- (sum(ttdf$adj.P.Val <= threshold)) / nrow(ttdf)
    }

  #Function def
    runIHWon1DF <- function(ttdf, alpha, proportion, ...){
      #run ihw on one df
      #return an ihwResult object
      if (is.null(ttdf$P.Value) || is.null(ttdf$AveExpr)){
        stop ("Expected both P.Value and AveExpr columns in data.frame")
      }
      IHWresult <- ihw(ttdf$P.Value,
                         covariates = ttdf$AveExpr,
                         alpha = alpha, ...)
    }

  #run IHW on each dataframe, collect the result objects in a list which
  #is added to the contrastList object.
  proportion <- sapply(contrastList$TopTableList, getProportion, threshold=FDRthreshold)
  contrastList$ihwList <- list()
  contrastNames <- names(contrastList$TopTableList)
  for (i in 1:length(contrastList$TopTableList)){
    message(paste("Running ihw on contrast", contrastNames[i], sep=" "))

    ihwResult = runIHWon1DF(contrastList$TopTableList[[i]],
                            alpha=alpha,
                            proportion=proportion[i], ...=...)
    #capture the ihwResult object
    contrastList$ihwList[[i]] <- ihwResult
    #cbind adj_pvalue, weight and weighted_pvalue columns to topTable DF
    contrastList$TopTableList[[i]] <- cbind (contrastList$TopTableList[[i]],
                          as.data.frame(ihwResult)[,2:4])
    #prefix the colnames of those three columns with "ihw."
    cnames = colnames(contrastList$TopTableList[[i]])
    numcol = length (cnames)
    cnames[(numcol-2):numcol] <- paste ("ihw.", cnames[(numcol-2):numcol], sep="")
    colnames(contrastList$TopTableList[[i]]) <- cnames
  }

  contrastList$workflowRecord$runIHW <- paste (date(), " : v", packageVersion("DGE.Tools"), sep="")

  #add optional facet plots of weights vs. Intensity colorby bin
  return(contrastList)

}
