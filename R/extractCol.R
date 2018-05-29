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
#' Note: all dataframes must have the same rownames/order
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords topTable
#'
#' @param dflist A list of data.frames which all have the same colnames and same row counts.
#' The dataframes in the list should have rownames (geneIDs).
#' @param ColName The name of the data column to extract to a matrix
#'
#' @return A matrix containing the extracted columns
#'
#' @examples
#' MyPvalues  = ExtractCol (MySLOA$TopTableList, "P.Value")
#'
#' @import magrittr
#'
#' @export
extractCol <- function(dflist, ColName){
  MyMatrix = lapply (dflist, `[[`, ColName) %>% do.call(what=cbind)
  #get gene ids from first df
  rownames(MyMatrix) = rownames(dflist[[1]])
  #transfer the contrast names
  colnames(MyMatrix) = names(dflist)
  return(MyMatrix)
}


#' @import magrittr dplyr tibble
#'
#' @export
extractCol2 <- function(dflist, ColName){
  #support combining topTable data from different DGEobjs
  for (i in 1:length(dflist)){

    newdat <- dflist[[i]] %>%
      rownames_to_column(var="rowid") %>%
      select(rowid, ColName)

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

