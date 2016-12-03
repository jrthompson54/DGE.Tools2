### Function extractCol ###
#' Function  extractCol
#'
#' Take a list of dataframes where each dataframe has the same
#' column names (e.g. a list of topTable dataframes). Extract
#' the named column from each dataframe and return a matrix.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords topTable
#'
#' @param dflist A list of data.frames which all have the same colnames.
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
