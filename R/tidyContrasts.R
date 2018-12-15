### Function tidyContrasts ###
#' Function  tidyContrasts
#'
#' Takes a DGEobj or contrast list as input and merges them into one dataframe.
#' A contrast list is simply a list of topTable contrast dataframes. For
#' example, DGEobj::getType(mydgeobj, "topTable") would retrieve a list of
#' contrasts from a DGEobj.  The contrast list must be a named list as the
#' contrast names are used during the merge operation.  Output can be specified
#' as "tall" (synonymous with "tidy") or "wide".
#'
#' The contrast names will be used as a new column in the tidy output format.
#' The input may or may not have rownames. If rownameColumn does not exist as a colname in
#' the dataframes, it is created from the rownames.  The output will always
#' contain a rownames column instead of actual rownames.
#'
#' The output will be in tall "tidy" format.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq; contrasts; merge
#'
#' @param x A DGEobj or named list of contrast dataframes
#'   (required).
#' @param rownameColumn Name of the rowname column. If a column by this
#'   name does not exist, it is created from the rownames property (rownames_to_column(var=rownameColumn))
#' @param columns A character vector of columns to include in the output (default = colnames of the first contrast)
#'
#' @return A DF with merged contrast data.
#'
#' @examples
#'
#'   MyMergedResults  = tidyContrasts (myDgeObj)
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr select bind_rows
#' @importFrom tidyr gather
#' @importFrom DGEobj getType
#'
#' @export
tidyContrasts <- function(x, columns, rownameColumn="rownames"){

  assertthat::assert_that(class(x)[[1]] %in% c("DGEobj", "list"))

  if (class(x)[[1]] == "DGEobj"){
    dgeObj <- x
    x <- DGEobj::getType(dgeObj, "topTable")
    if (length(x) == 0)
      stop("No topTable dataframes found in DGEobj")
  }  #x is now a contrastlist

  #make sure list contains only dataframes
  if (all(sapply(x, class) == "data.frame") == FALSE)
    stop ("Input list must contain only dataframes\n")

  #make sure each df has a name.
  minNameLen <- min(sapply(names(x), nchar))
  if (minNameLen == 0) {
    stop("At least one of the dataframes in the input list has no name.\n")
  }

  #Set default columns
  if (missing(columns))
    columns <- colnames(x[[1]])

  #find the common set of columns present in all dataframes.
  commonColumns <- colnames(x[[1]])
  for (i in 2:length(x))
    commonColumns <- intersect(commonColumns, colnames(x[[i]]))

 #make sure user-requested columns are present
  if (!all(columns %in% commonColumns))
    warning("Some requested columns are not present in all dataframes.")
  commonColumns <- intersect(commonColumns, columns)

  #Does the rownameColumn exist in df1?
  if (!rownameColumn %in% colnames(x[[1]])) {
    #move rownames to rownameColumn
    x <- lapply(x, rownames_to_column, var=rownameColumn)
    commonColumns <- c(rownameColumn, commonColumns)
  }

  #reduce all dataframes to the selected columns
  x <- lapply(x, select, commonColumns)

  #add a contrast name column to each DF
  for (name in names(x)){
    x[[name]]["contrast"] <- name
  }

  #Now merge the dataframes vertically
  x <- dplyr::bind_rows(x)

  return(x)
}



