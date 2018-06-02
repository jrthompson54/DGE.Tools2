### Function extractCol ###
#' Function  extractCol
#'
#' Take a named list of dataframes where each dataframe has the same
#' column names (e.g. a list of topTable dataframes). Extract
#' the named column from each dataframe and return a matrix.
#'
#' The common use case for this is to provide a list of topTable
#' data frames and extract one column from each file to create
#' a matrix of LogRatios or Pvalues.
#'
#' Technically, this should work as long as the requested colName is present
#' in each dataframe.  The default robust = TRUE should be used unless you
#' are absolutely certain each dataframe in the input list has the same row count
#' and row order.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords topTable
#'
#' @param dflist A list of data.frames which all have the same colnames and same row counts.
#' The dataframes in the list should have rownames (geneIDs).
#' @param colName The name of the data column to extract to a matrix
#' @param robust Default = TRUE;  TRUE forces use of a joins to merge columns
#'   which is more reliable, allows you to combine contrasts from different
#'   projects but may nor return items in the same row order as the source
#'   table.  Setting to false invokes a cbind approach that requires all
#'   dataframes to have the same row count and row order but preserves the
#'   original row order
#'
#' @return A matrix containing the extracted columns
#'
#' @examples
#'
#'   MyPvalues  = ExtractCol (TopTableList, colName="P.Value")
#'
#' @import magrittr
#'
#' @export
extractCol <- function(dflist, colName, robust="TRUE"){

  ifelse(robust,
         .extracCol2(dflist, colName),
         .extracCol1(dflist, colName)
  )
}

### Function extractCol1 ###
#' Function  extractCol1
#'
#' Take a named list of dataframes where each dataframe has the same
#' column names (e.g. a list of topTable dataframes). Extract
#' the named column from each dataframe and return a matrix.
#'
#' The common use case for this is to provide a list of topTable
#' data frames and extract one column from each file to create
#' a matrix of LogRatios or Pvalues.
#'
#' Note: all dataframes must have the same rownames/order
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords topTable
#'
#' @param dflist A list of data.frames which all have the same colnames and same row counts.
#' The dataframes in the list should have rownames (geneIDs).
#' @param colName The name of the data column to extract to a matrix
#'
#' @return A matrix containing the extracted columns
#'
#' @examples
#' MyPvalues  = ExtractCol1 (TopTableList, colName="P.Value")
#'
#' @import magrittr
#'
.extractCol1 <- function(dflist, colName){
  MyMatrix = lapply (dflist, `[[`, colName) %>% do.call(what=cbind)
  #get gene ids from first df
  rownames(MyMatrix) = rownames(dflist[[1]])
  #transfer the contrast names
  colnames(MyMatrix) = names(dflist)
  return(MyMatrix)
}

### Function extractCol2 ###
#' Function  extractCol2
#'
#' Take a named list of dataframes where each dataframe has the same
#' column names (e.g. a list of topTable dataframes). Extract
#' the named column from each dataframe and return a matrix.
#'
#' The common use case for this is to provide a list of topTable
#' data frames and extract one column from each file to create
#' a matrix of LogRatios or Pvalues.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords topTable
#'
#' @param dflist A list of data.frames which all have the same colnames and same row counts.
#' The dataframes in the list should have rownames (geneIDs).
#' @param colName The name of the data column to extract to a matrix
#'
#' @return A matrix containing the extracted columns
#'
#' @examples
#' MyPvalues  = ExtractCol TopTableList, "P.Value")
#'
#' @import magrittr dplyr tibble
.extractCol2 <- function(dflist, colName){
  #support combining topTable data from different DGEobjs
  for (i in 1:length(dflist)){

    newdat <- dflist[[i]] %>%
      rownames_to_column(var="rowid") %>%
      select(rowid, colName)

    if (i == 1){
      dat <- newdat
    } else {
      dat %<>% full_join(newdat, by="rowid")
    }
  }
  dat %<>% column_to_rownames(var="rowid")
  colnames(dat) <- names(dflist)
  return(dat)
}

