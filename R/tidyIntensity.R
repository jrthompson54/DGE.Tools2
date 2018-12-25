### Function tidyIntensity ###
#' Function  tidyIntensity
#'
#' Takes dataframe or matrix of intensity values with geneid rownames and samplename columnnames and a grouping vector.
#' Returns a tidy dataframe.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq; intensity; tidy
#'
#' @param x A dataframe or matrix with rownames and colnames (required).
#' @param rowidColname Column name for the rowID (Usually a geneID column).
#'   If it doesn't exist, it will be created from rownames.
#' @param valueColname Column name for the value column
#' @param group A grouping variable to define which columns go together (vector
#'   with length = ncol(x); coercable to character)
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
}


