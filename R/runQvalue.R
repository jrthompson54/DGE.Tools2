### Function runQvalue ###
#' Function  runQvalue
#'
#' Takes an list of contrasts (e.g. topTable output or other dataframes that contain
#' a pvalue column).  Adds a qvalue and local FDR (lfdr) column to each dataframe.
#'
#' The qvalue package from John Storey at Princeton takes a list of pvalues and
#' calculates a qvalue and a Local FDR (LFDR). The qvalue is essentially a less
#' conservative FDR estimate compared to the default Benjamini-Hochberg FDR
#' produced by topTable analysis (i.e. will give you more differential genes at
#' the same nominal cutoff). The qvalue function also produces a Local FDR
#' (lfdr) column which answers a slightly different and possibly more relevant
#' question.  The BH FDR (adj.P.Val in topTable data.frames) and qvalue tell you
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
#' @param contrastList An list of dataframes with a pvalue column (all tables
#'   must use the same colname for the pval column)
#' @param pvalField Optional. (Default = "P.Value")  Not needed if you're working
#'   with topTable output.  Use to define the colname of the pvalue field in
#'   each dataframe.
#' @param ... Optional arguments passed to the qvalue function (See ?qvalue)
#'
#' @return The input contrastList now containing qvalue and lfdr columns
#'   in each dataframes.
#'
#' @examples
#' MyContrastList <- runQvalue(MyContrastList)
#' #The magrittr way:
#' MyContrastList %<>% runQvalue
#'
#' @import qvalue magrittr
#'
#' @export
runQvalue <- function(contrastList, pvalField="P.Value", ...){
### Add Qvalues to each topTable dataframe in contrastList ###

  assertthat(class(contrastList)[[1]] == "list")

  contrastNames = names(contrastList)

  for (i in 1:length(contrastList)) {
    assertthat(exists(pvalField, contrastList[[i]]))
    p = contrastList[[i]][, pvalField]
    q = qvalue(p, ...)
    #add the qvalue and lfdr columns to the topTable df
    contrastList[[i]]$Qvalue = q$qvaluesq
    contrastList[[i]]$qvalue.lfdr = q$lfdr
    #add documentation
    attr(contrastList[[i]], "qvalue") = TRUE
  }

  return(contrastList)
}



