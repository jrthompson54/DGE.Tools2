### Function tidyIntensity ###
#' Function  tidyIntensit
#'
#' Takes dataframe or matrix of intensity values with geneid rownames and samplename columnnames.
#'
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
#' @param x A dataframe or matrix with rownames (required).
#' @param rowIDcolumn Name of column for the rowID (Usually a geneID column).  If it doesn't exist, it will be created from rownames (default="rowname)
#' @param group A grouping variable to define which columns go together (vector with length = ncol(x); coercable to character)
#'
#' @return A tidy dataframe of intensity data.
#'
#' @examples
#'
#'   MyTidyIntensity  = tidyIntensity (MyLog2CPM)
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#'
#' @export
tidyIntensity <- function(x, rowIDcolumn="rownames", group){

  assertthat::assert_that(class(x)[[1]] %in% c("data.frame", "matrix"),
                          !missing(group),
                          length(group) == ncol(x))

  if (class(x)[[1]] == "matrix") x <- as.data.frame(x)
  sample <- colnames(x)

  #create a rownames column
  x <- tibble::rownames_to_column(x, var=rowIDcolumn)

  x <- gather(x, key="sample", value="log2cpm", -rowIDcolumn)

  #join the group info
  groupinfo <- data.frame(sample=sample, group=group, stringsAsFactors =FALSE)
  x <- dplyr::left_join(x, groupinfo, by="sample")

}



