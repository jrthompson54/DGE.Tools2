### Function runIHW ###
#' Function  runIHW
#'
#' This is a wrapper around the independent hypothesis weighting package that
#' takes a list of topTable dataframes and applies Independent Hypothesis
#' Weighting (IHW) to each topTable dataframe in the list.
#'
#' IHW is a method developed by N. Ignatiadis (http://dx.doi.org/10.1101/034330)
#' to weight FDR values based on a covariate (AveExpr in this case).
#'
#' The IHW FDR values are added as additional columns to the topTable dataframes.
#'
#' Note that it is impractical to run IHW on a list of genes less than ~5000
#' because operationally, it breaks the data into bins of 1500 genes for the
#' analysis if bins = 1, IHW converges on the BH FDR value and there's no point.
#' So we run IHW on the complete set of detected genes from topTable results
#' (not topTreat).
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2, png, bmp, tiff, jpeg, pdf
#'
#' @param contrastList A list of topTable dataframes.
#' @param alpha The alpha value is the FDR level you wish to interogate (range 0-1; default = 0.1)
#' @param ... other arguments are passed directly to the ihw function (see ?ihw).
#'
#' @return A list of lists.  The first element is the original contrastList with
#'   additional IHW columns added to each dataframe.   The TopTable dataframes
#'   will contain additional columns added by the ihw analysis and prefixed with
#'   "ihw." The second list element is the IHW result dataframe.
#'
#' @examples
#' IHWresults <- runIHW(MyContrastList)
#' MyContrastList <- IHWresults[[1]]
#' IHWdf <- IHWresults[[2]]
#'
#' @importFrom IHW ihw
#'
#' @export
 runIHW <- function(contrastList,
                    alpha=0.1,
                    FDRthreshold = 0.1,
                    ...){

  #Function def
    getProportion <- function(ttdf, threshold) {
      #get the proportion for one df
      bhfdrproportion <- (sum(ttdf$adj.P.Val <= threshold)) / nrow(ttdf)
    }

    #note proportion is passed below in runIHWon1DF but not used!?!?
    #Leave for now but consider trimming later.
  #Function def
    runIHWon1DF <- function(ttdf, alpha, proportion, ...){
      #run ihw on one df
      #return an ihwResult object
      if (is.null(ttdf$P.Value) || is.null(ttdf$AveExpr)){
        stop ("Expected both P.Value and AveExpr columns in data.frame")
      }
      IHWresult <- IHW::ihw(ttdf$P.Value,
                         covariates = ttdf$AveExpr,
                         alpha = alpha, ...)
    }

  #run IHW on each dataframe, collect the result objects in a list which
  #is added to the result object.
  proportion <- sapply(contrastList, getProportion, threshold=FDRthreshold)
  ihwList <- list()
  contrastNames <- names(contrastList)
  for (i in 1:length(contrastList)){
    # message(paste("Running ihw on contrast", contrastNames[i], sep=" "))

    ihwResult <- runIHWon1DF(contrastList[[i]],
                            alpha=alpha,
                            proportion=proportion[i], ...)
    #capture the ihwResult object
    ihwList[[i]] <- ihwResult
    #cbind adj_pvalue, weight and weighted_pvalue columns to topTable DF
    contrastList[[i]] <- cbind (contrastList[[i]],
                          as.data.frame(ihwResult)[,2:4])
    #prefix the colnames of those three columns with "ihw."
    cnames <- colnames(contrastList[[i]])
    numcol <- length (cnames)
    cnames[(numcol-2):numcol] <- paste ("ihw.", cnames[(numcol-2):numcol], sep="")
    colnames(contrastList[[i]]) <- cnames

    #add documentation
    attr(contrastList[[i]], "ihw") = TRUE
  }

  result <- list(contrasts=contrastList, ihwObj=ihwList)
  return(result)

}
