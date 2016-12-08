### DGEresult object definition

### # Constructor function for the class
DGEresult <- function(item, itemName, itemType) {
    #initialize an empty DGEresult and optionally add the first item
    dgeresult <- list()
    dgeresult$data <- list()
    dgeresult$type <- list()
    class(dgeresult) <- "DGEresult"
    if (!missing(item) & !missing(itemName) & !missing(itemType))
        dgeresult <- addItem.DGEresult(dgeresult, item, itemName, itemType)
    return (dgeresult)
}

### .checkDimension
.checkDimension <- function(x, ...) UseMethod(".checkDimension")
.checkDimension.default <- function(dgeResult, ...)
{
    warning(paste(".checkDimension does not know how to handle object of class ", 
                  class(dgeResult), 
                  "and can only be used on class DGEresult"))
}
.dimensionMatch.DGEresult <- function(dgeResult, item, itemType){
    #check the dimensions of item against the DGEresult object dim 
    #set the DGEresult dimension if not set (i.e. case of empty DGEresult)
    #stop with an error if item is mismatched with DGEresult dim
    result <- FALSE
    switch(itemType,
           "row" = {
               if (dim(dgeResult)[1] >0 & dim(dgeResult)[1] != nrow(item))
                   stop('New row object does not match row dimension of DGEresult object')
           },
           
           "col" = {
               if (dim(dgeResult)[2] >0 & dim(dgeResult)[2] != nrow(item))
                   stop('Rows in new col object does not match col dimension of DGEresult object')
           },
           
           "assay" = {
               if (dim(dgeResult)[1] >0 & dim(dgeResult)[1] != nrow(item))
                   stop('New assay object does not match row dimension of DGEresult object')
               if (dim(dgeResult)[2] >0 & dim(dgeResult)[2] != nrow(item))
                   stop('New assay object does not match col dimension of DGEresult object')
           },
           
           "contrastTop" = {
               if (dim(dgeResult)[1] >0 & dim(dgeResult)[1] != nrow(item))
                   stop('New contrastTop object does not match row dimension of DGEresult object')
           }
           
    ) #end switch
    result <- TRUE
    return(result)
}

### addItem
# addItem <- function(dgeResult, item, itemName, itemType,
#                     overwrite=FALSE)
# {
#     #print("Calling the base addItem function")
#     UseMethod("addItem", dgeResult)
# }
addItem <- function(x, ...) UseMethod("addItem")

# addItem.default <- function(dgeResult, item, itemName, itemType,
#                             overwrite=FALSE)
# 
addItem.default <- function(dgeResult, ...)
{
    warning(paste("addItem does not know how to handle object of class ", 
                  class(dgeResult), 
                  "and can only be used on class DGEresult"))
}
addItem.DGEresult <- function(dgeResult, item, itemName, itemType,
                              overwrite=FALSE){
    #enforce itemType
    if (!itemType %in% c("row", "col", "assay", "contrastTop", "contrastTreat"))
        stop("itemType must be one of: row, col, assay, contrastTop, contrastTreat")
    
    #add item and assign it's type
    #refuse to add if it exists already unless overwrite = T
    if (overwrite==FALSE & with(dgeResult$data, exists(itemName)))
        stop('itemName already exists in DGEresult!')
    #confirm dimensions consistent before adding
    else if (.dimensionMatch(dgeResult, item, itemType) == TRUE){
        dgeResult$data[[itemName]] <- item
        dgeResult$type[[itemName]] <- itemType
    }
    return(dgeResult)  
}

### rmItem
rmItem <- function(x, ...) UseMethod("rmItem")
rmItem.default <- function(dgeResult, ...)
{
    warning(paste("rmItem does not know how to handle object of class ", 
                  class(dgeResult), 
                  "and can only be used on class DGEresult"))
}
rmItem.DGEresult <- function(dgeResult, itemName){
    #remove the named item
    dgeResult$data[[itemName]] <- NULL
    dgeResult$type[[itemName]] <- NULL
}

### getItem
getItem <- function(x, ...) UseMethod("dim")
getItem.default <- function(dgeResult, ...)
{
    warning(paste("getItem does not know how to handle object of class ", 
                  class(dgeResult), 
                  "and can only be used on class DGEresult"))
}
getItem.DGEresult <- function(dgeResult, itemName){
    return(dgeResult$data[[itemName]])
}

### dim
dim <- function(x, ...) UseMethod("dim")
dim.default <- function(dgeResult, ...)
{
    warning(paste("dim does not know how to handle object of class ", 
                  class(dgeResult), 
                  "and can only be used on class DGEresult"))
}
dim.DGEresult <- function(dgeResult){
    
    dimension <- c(0,0)
    rowcount <- 0
    colcount <- 0
    #look for first row or assay or contrastTop object and get nrow
    idx <- dgeResult[["type"]] %in% c("row", "assay", "contrastTop")
    if (sum(idx) > 0) {
        firstItemName <- names(dgeResult$type[idx])[[1]]
        rowcount <- nrow(dgeResult$data[[firstItemName]])
    }
    
    #look for first col or assay element
    idx <- dgeResult[["type"]] %in% c("col", "assay")
    if (sum(idx) > 0) {
        firstItemName <- names(dgeResult$type[idx])[[1]]
        firstItemType <- dgeResult$type[[firstItemName]]
        firstItem <- dgeResult$data[[firstItemName]]
        if (firstItemType == "col")
            colcount <- nrow(firstItem)
        else colcount <- ncol(firstItem)
    }
    dimension <- c(rowcount, colcount)
    return(dimension)
}

