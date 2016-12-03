### Function runQvalue ###
#' Function  runQvalue
#'
#' Takes an SLOA object from runContrasts and add a qvalue and local FDR (lfdr)
#' column to each topTable result in the SLOA$TopTableList. Returns the modified
#' SLOA object.
#'
#' The qvalue package from John Storey at Princeton takes a list of pvalues and
#' calculates a qvalue and a Local FDR (LFDR). The qvalue is essentially a less
#' conservative FDR estimate compared to the default Benjamini-Hochberg FDR
#' produced by topTable analysis (i.e. will give you more differential genes at
#' the same nominal cutoff). The qvalue function also produces a Local FDR
#' (lfdr) column which answers a slightly different and possibly more relevant
#' question.  The BHFDR (adj.P.Val in topTable data.frames) and qvalue tell you
#' what the false discovery rate is for a list of genes at a given threshold.
#' The local FDR attempts to answer the question: what is the probability that
#' this particular gene is a false discovery?
#'
#' See \url{http://genomine.org/papers/Storey_FDR_2011.pdf} for a brief introduction
#' to FDRs and q-values.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords qvalue; lfdr; contrasts
#'
#' @param contrastList An SLOA object resulting from runContrasts
#' @param ... Optional arguments passed to the qvalue function (See ?qvalue)
#'
#' @return The modified input contrastList now containing qvalue and lfdr data
#'   in the topTable dataframes.
#'
#' @examples
#' MyContrastList <- runQvalue(MyContrastList)
#' #The magrittr way:
#' MyContrastList %<>% runQvalue
#'
#' @import qvalue magrittr
#'
#' @export
runQvalue <- function(contrastList, ...){
### Add Qvalues to each topTable dataframe in contrastList$TopTableList ###

  FVersion <- "07Jan2016"

  Qlist <- list()
  contrastNames = names(contrastList$TopTableList)

  for (i in 1:length(contrastList$TopTableList)) {
    q = qvalue(contrastList$TopTableList[[i]]$P.Value, ...)
    #add the qvalue and lfdr columns to the topTable df
    contrastList$TopTableList[[i]]$Qvalue = q$qvalues
    contrastList$TopTableList[[i]]$qvalue.lfdr = q$lfdr
    #capture the qvalue object
    Qlist[[i]] <- q
  }

  #name the qvalue objects
  names(Qlist) <- contrastNames
  contrastList$Qlist <- Qlist
  contrastList$workflowRecord$runQvalue <- paste (date(), " : v", packageVersion("DGE.Tools"), sep="")
  contrastList$workflowRecord$Date = date()

  return(contrastList)
}



