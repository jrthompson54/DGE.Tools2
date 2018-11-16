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
#' @return A dataframe containing the extracted columns
#'
#' @examples
#'
#'   MyPvalues  = ExtractCol (TopTableList, colName="P.Value")
#'
#' @export
extractCol <- function(dflist, colName, robust="TRUE"){

  ifelse(robust,
         return(.extractCol2(dflist, colName)),
         return(.extractCol1(dflist, colName))
  )
}

### Function .extractCol1 ###
#' Function  .extractCol1
#'
#' Take a named list of dataframes where each dataframe has the same
#' column names (typically a list of topTable dataframes). Extract
#' the named column from each dataframe and return a matrix.
#'
#' The common use case for this is to provide a list of topTable
#' data frames and extract one column from each file to create
#' a matrix of LogRatios or Pvalues.
#'
#' NOTE: all dataframes must have the same rownames and roworder
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords topTable
#'
#' @param dflist A list of data.frames which all have the same colnames and same row counts.
#' The dataframes in the list should have rownames (geneIDs).
#' @param colName The name of the data column to extract to a matrix
#'
#' @return A dataframe containing the extracted columns
#'
#' @examples
#' MyPvalues  = ExtractCol1 (TopTableList, colName="P.Value")
#'
#' @import magrittr
#' @importFrom assertthat assert_that
#'
.extractCol1 <- function(dflist, colName){
  assertthat::assert_that(class(dflist)[[1]] == "list",
                          class(colName)[[1]] == "character",
                          !is.null(names(dflist)))

  MyMatrix = lapply (dflist, `[[`, colName) %>% do.call(what=cbind)
  #get gene ids from first df
  rownames(MyMatrix) <- rownames(dflist[[1]])
  #transfer the contrast names
  colnames(MyMatrix) <- names(dflist)
  return(as.data.frame(MyMatrix))
}

### Function .extractCol2 ###
#' Function  .extractCol2
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
#' @param dflist A named list of dataframes which all have some common colnames and
#'   some overlap on rownames (genes). and same row counts. The dataframes in
#'   the list should have rownames (geneIDs).
#' @param colName The name of the data column to extract to a matrix.  This
#'   column must be present in all dataframes on the list.
#'
#' @return A dataframe containing the extracted columns
#'
#' @examples
#'
#'     MyPvalues  = .extractCol2(TopTableList, "P.Value")
#'
#' @import magrittr
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr select full_join
#' @importFrom assertthat assert_that
#'
.extractCol2 <- function(dflist, colName){
  #support combining topTable data from different DGEobjs

  assertthat::assert_that(class(dflist)[[1]] == "list",
                          class(colName)[[1]] == "character",
                          !is.null(names(dflist)))

  for (i in 1:length(dflist)){

    newdat <- dflist[[i]] %>%
      tibble::rownames_to_column(var="rowid") %>%
      dplyr::select(rowid, colName)

    if (i == 1){
      dat <- newdat
    } else {
      dat %<>% dplyr::full_join(newdat, by="rowid")
    }
  }
  dat %<>% tibble::column_to_rownames(var="rowid")
  colnames(dat) <- names(dflist)
  return(dat)
}

