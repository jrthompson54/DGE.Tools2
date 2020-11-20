### Function df_to_excel ###
#' Function  df_to_excel
#'
#' Take a dataframe with rownames and output to an Excel spreadsheet.
#' Freezes the row and colnames and sets columns into Excel filter mode.
#'
#' @author John Thompson, \email{jrt@@thompsonclan.org}
#' @keywords dataframe, Excel
#'
#' @param df A dataframe with rownames and colnames (Required).
#' @param filename Path/Filename for output (Required).
#' @param sheetname Name for the Excel sheet (Default = "data")
#' @param widthFirstCol width of 1st column (rownames) (Default=40)
#' @param widthOtherCol width of 1st column (rownames) (Default=10)
#' @param baseFontSize base font for Excel file (Default=10)
#' @param useRownames Creates a spreadsheet without rownames (default=TRUE).  Set to FALSE
#'   when there is no column suitably unique to use for rownames
#'
#' @return Writes an Excel file to the specified file.
#'
#' @examples
#' #merge a set of topTable contrasts into one data.frame
#' MyContrastTable  = topTable.merge (topTablelist)
#' df_to_excel(MyContrastTable, "myExcelFile.xlsx")
#'
#' @import magrittr
#' @importFrom stringr str_c
#' @importFrom assertthat assert_that
#' @importFrom openxlsx createWorkbook modifyBaseFont addWorksheet freezePane writeDataTable setColWidths saveWorkbook
#'
#' @export
df_to_excel <- function(df, filename,
                        sheetname="data",
                        widthFirstCol = 40,
                        widthOtherCol = 10,
                        baseFontSize = 10,
                        useRownames = TRUE
                        ){

  assertthat::assert_that(!missing(df),
                          "data.frame" %in% class(df),
                          !missing(filename),
                          class(filename)[[1]] == "character"
                          )

  numcol <- ncol(df)
  wb <- openxlsx::createWorkbook()
  options("openxlsx.borderColour" = "#4F80BD")
  options("openxlsx.borderStyle" = "thin")
  openxlsx::modifyBaseFont(wb, fontSize = baseFontSize, fontName = "Arial Narrow")

  openxlsx::addWorksheet(wb, sheetName = sheetname, gridLines = TRUE)
  openxlsx::freezePane(wb, sheet = sheetname, firstRow = TRUE, firstCol = TRUE)
  openxlsx::writeDataTable(wb, sheet = sheetname, x = df, colNames = TRUE, rowNames = useRownames, tableStyle = "TableStyleLight9")
  openxlsx::setColWidths(wb, sheet = sheetname, cols = 1, widths = widthFirstCol)
  openxlsx::setColWidths(wb, sheet = sheetname, cols = 2:numcol, widths = widthOtherCol)

  openxlsx::saveWorkbook(wb, filename, overwrite = TRUE)

}


### Function dflist_to_excel ###
#' Function  dflist_to_excel
#'
#' Take a named list of dataframes with rownames and output to an Excel spreadsheet.
#' Freezes the row and colnames and sets columns into Excel filter mode.
#'
#' @author John Thompson, \email{jrt@@thompsonclan.org}
#' @keywords dataframe, Excel
#'
#' @param dflist A named list of dataframes with rownames and colnames (Required).
#' @param filename Path/Filename for output (Required).
#' @param widthFirstCol width of 1st column (rownames) (Default=40)
#' @param widthOtherCol width of 1st column (rownames) (Default=10)
#' @param baseFontSize base font for Excel file (Default=10)
#' @param useRownames default = TRUE; Set to false when no appropriate unique column present.
#'
#' @return Writes an Excel file to the specified file.
#'
#' @examples
#' MyContrastTable  = topTable.merge (topTablelist)
#'
#' @import magrittr
#' @importFrom stringr str_c
#' @importFrom assertthat assert_that
#' @importFrom openxlsx createWorkbook modifyBaseFont addWorksheet freezePane writeDataTable setColWidths saveWorkbook
#'
#' @export
dflist_to_excel <- function(dflist, filename,
                        widthFirstCol = 40,
                        widthOtherCol = 10,
                        baseFontSize = 10,
                        useRownames = TRUE
){

  .addsheet <- function(wb, df, sheetname, widthFirstCol, widthOtherCol, useRownames=useRownames){
    openxlsx::addWorksheet(wb, sheetName = sheetname, gridLines = TRUE)
    openxlsx::freezePane(wb, sheet = sheetname, firstRow = TRUE, firstCol = TRUE)
    openxlsx::writeDataTable(wb, sheet = sheetname, x = df, colNames = TRUE, rowNames = useRownames, tableStyle = "TableStyleLight9")
    openxlsx::setColWidths(wb, sheet = sheetname, cols = 1, widths = widthFirstCol)
    openxlsx::setColWidths(wb, sheet = sheetname, cols = 2:ncol(df), widths = widthOtherCol)
    return(wb)
  }

  assertthat::assert_that(!missing(dflist),
                          class(dflist)[[1]] == "list",
                          !is.null(names(dflist)),
                          !missing(filename),
                          class(filename)[[1]] == "character"
  )

  numcol <- ncol(df)
  wb <- openxlsx::createWorkbook()
  options("openxlsx.borderColour" = "#4F80BD")
  options("openxlsx.borderStyle" = "thin")
  openxlsx::modifyBaseFont(wb, fontSize = baseFontSize, fontName = "Arial Narrow")

  for (name in names(dflist)){
    wb <- .addsheet(wb, dflist[[name]], name, widthFirstCol, widthOtherCol)
  }

  openxlsx::saveWorkbook(wb, filename, overwrite = TRUE)

}

