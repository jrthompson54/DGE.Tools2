### Function tidyIntensity ###
#' Function  tidyIntensity
#'
#' Takes dataframe or matrix of intensity values with geneid rownames and samplename columnnames and a grouping vector.
#' Returns a tidy dataframe.
#'
#' This function takes a genes x samples datafile of e.g. log2cpm.  It uses tidyr::gather to produce
#' a "tidy" datafile.
#'
#' So the rowidColname defines the name you want the gene id column to have and will create the column from
#' the rownames which are typically ensembl gene IDs.  If you add a gene symbol column to your intensity
#' dataframe, you can specify the name of that column with the rowidColname argument and it will be used instead
#' of the rownames.
#'
#'
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq; intensity; tidy
#'
#' @param x A dataframe or matrix with rownames and colnames (required)
#' @param rowidColname Column name for the rowID (Usually a geneID column)
#'   If it doesn't exist, it will be created from rownames  (required)
#' @param keyColname Defines the column name that will hold the sample names
#'   captured from the colnames of the input dataframe  (required)
#' @param valueColname Defines the column name for the value column (intensity values)  (required)
#' @param group A grouping variable to define which columns go together (vector
#'   with length = ncol(x) of the input dataframe).  This is typically the name
#'   of a column from the design table  (required)
#'
#' @return A tidy dataframe of intensity data.
#'
#' @examples
#'
#'   MyTidyIntensity  = tidyIntensity (MyLog2CPM,
#'                                     rowidColname="GeneID",
#'                                     keyColname="Sample"
#'                                     valueColname="Log2CPM",
#'                                     group=dgeObj$design$ReplicateGroup)
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#'
#' @export
tidyIntensity <- function(x, rowidColname, keyColname, valueColname, group){

  assertthat::assert_that(class(x)[[1]] %in% c("data.frame", "matrix"),
                          !missing(rowidColname),
                          !missing(keyColname),
                          !missing(valueColname),
                          !missing(group),
                          length(group) == ncol(x))

  if (class(x)[[1]] == "matrix") x <- as.data.frame(x)

  samples <- colnames(x)

  #create a rownames column
  xcol <- ncol(x)+1
  x <- tibble::rownames_to_column(x, var=rowidColname)
  x <- tidyr::gather(x, key=!!keyColname, value=!!valueColname, 2:xcol)

  #join the group info
  # groupinfo <- tibble(!!keyColname:=sample, group=group)
  # x <- dplyr::left_join(x, groupinfo, by=keyColname)
  groupinfo <- tibble(keycol=samples, group=group)
  x <- dplyr::left_join(x, groupinfo, by=rlang::set_names(nm = keyColname, "keycol"))
  return(x)
}


