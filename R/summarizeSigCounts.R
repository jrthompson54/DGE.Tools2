### Function  summarizeSigCounts ###
#' Function  summarizeSigCounts
#'
#' Takes a contrast list produced by runContrasts. Also, specify columns to
#' summarize and thresholds for each column.  Optionally, provide a fold change
#' threshold.  Queries the topTable results and returns a dataframe with the
#' summary results; gene counts that meet the specified conditions.
#'
#' Specified column names that don't exist will be ignored. So normally, the
#' defaults cover all the pvalue and FDR related columns and you'll
#' only have to add a fcThreshold if you want one or modify the pvalue/fdr
#' thresholds if you want different thresholds from the defaults.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords DGE.Tools, signature summary, gene counts
#'
#' @param ContrastList A list of topTable dataframes
#' @param columns Vector of column names to summarize from topTable dataframes
#'   Default = c("P.Value", "adj.P.Val", "Qvalue", "qvalue.lfdr", "ihw.adj_pvalue")
#' @param sigThresholds Thresholds to use for each column specified in columns
#'   Must be same length at columns argument.
#'   Default = c(0.01, 0.1, 0.1, 0.1, 0.1)
#' @param fcThreshold vector of thresholds for specified fields
#'
#' @return data.frame with one summary row per contrast
#'
#' @import magrittr
#'
#' @examples
#' 
#' #get a contrast list from a dgeObj
#' MyContrastList <- getType(dgeObj, "topTable")
#'
#' #all default thresholds, no fold change threshold
#' MySigSummary <- summarizeSigCounts (MyContrastList)
#'
#' #all defaults with a foldchange threshold
#' MySigSummary <- summarizeSigCounts (MyContrastList, fcThreshold=1.5)
#'
#' #change the pvalue and fdr thresholds
#' MySigSummary <- summarizeSigCounts (MyContrastList,
#'     sigThresholds = c(0.05, 0.2, 0.2, 0.2, 0.2))
#'
#' @export
summarizeSigCounts <- function(ContrastList,
  columns = c("P.Value", "adj.P.Val", "Qvalue", "qvalue.lfdr", "ihw.adj_pvalue"),
  sigThresholds = c(0.01, 0.1, 0.1, 0.1, 0.1),
  fcThreshold = 0
  ){

  #Functions
  getSigCounts <- function(df,
                           columns = c("P.Value", "adj.P.Val", "Qvalue", "qvalue.lfdr", "ihw.adj_pvalue"),
                           thresholds = c(0.01, 0.1, 0.1, 0.1, 0.1),
                           fcThreshold=0){
    #pass a topTable df, vector of field names, a threshold and optionally a fcThreshold)
    #Then get the sigcounts
    #return a vector of named values
    CountSig <- function(df, column, threshold, fcThreshold=0){
      # supply a field and a threshold.
      # return the number of samples <= threshold
      # If fcThreshold >0, also filter on logFC field (convert the fcThreshold to
      # log2 and filter the abs of logFC)

      idx <- df[column] <= threshold
      if (fcThreshold>0) {
        fcidx <- abs(df$logFC) >= log2(fcThreshold)
        idx <- idx & fcidx
      }
      return(sum(idx))
    }  #End F. CountSig

    #body getSigCounts
    counts <- list()

    for (i in 1:length(columns)){
      counts[columns[i]] <- CountSig(df, columns[i], thresholds[i], fcThreshold)
    }
    return(unlist(counts))
  } #end F. getSigCounts

  #Body summarizeSigCounts

  if (!length(columns) == length(sigThresholds)){
    stop("sigThresholds argument should be same length as columns argument!")
  }

  #reduce tableFields to only ones that exist
  columns <- columns[columns %in% colnames(ContrastList[[1]])]
  sigThresholds <- sigThresholds[columns %in% colnames(ContrastList[[1]])]

  #collect the rows in a list
  myrows <-list()
  for (i in 1:length(ContrastList)){
    myrows[[i]] <- getSigCounts(ContrastList[[i]], columns, sigThresholds, fcThreshold)
  }

 #put rows into a matrix
 DF <- do.call(rbind, myrows)

 rownames(DF) <- names(ContrastList)
 colnames(DF) <- columns

 return(DF)
}




