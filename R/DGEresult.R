### DGEresult object definition


### Define an extensible datatype
### There are 4 imutable base type: row, col, assay, meta
### Each type must be one of these four basetypes
### Additional types can be added to .type as long as they are
### assigned to one of the 4 basetypes.
.basetype <- c("row", "col", "assay", "meta")
#the value of .type is a basetype
.type <- list(row="row", col="col", assay="assay", meta="meta",
               contrastTop = "row", contrastTreat="meta")
# x = getwd()
# setwd ("~/R/lib/pkgsrc/DGE.Tools2/")
# save(.basetype, file="./data/basetype.rda")
# save(.type, file="./data/type.rda")
# setwd(x)



### # Constructor function for the class
DGEresult <- function(item, itemName, itemType, workflow) {
    #initialize an empty DGEresult and optionally add the first item
    #Workflow is optional and should be a list with 4 specific elements:
    #   workflow$itemName <- list(itemName=itemName,        
    #                           type=type, 
    #                           Function=callingFunctionName,
    #                           Date=date()
    #                           Function=funName
    #                           args=FunArgs)

    dgeresult <- list()
    dgeresult$data <- list()
    dgeresult$type <- list()
    dgeresult$workflow <- list()
    dgeresult$workflow[["DGEresult"]] <- list(
        itemName = "DGEresult",
        type = "DGEresult",
        date = date(),
        funArgs = match.call()
    )
    class(dgeresult) <- "DGEresult"
    if (!missing(item) & !missing(itemName) & !missing(itemType))
        dgeresult <- addItem.DGEresult(dgeresult, item, itemName, itemType)
    if (!missing(workflow)){
    }
    return (dgeresult)
}

### .checkDimension
.dimensionMatch <- function(x, ...) UseMethod(".dimensionMatch")
.dimensionMatch.default <- function(dgeResult, ...) {
    warning(paste(".dimensionMatch does not know how to handle object of class ", 
                  class(dgeResult), 
                  "and can only be used on class DGEresult"))
}
.dimensionMatch.DGEresult <- function(dgeResult, item, itemType){
    #check the dimensions of item against the DGEresult object dim 
    #set the DGEresult dimension if not set (i.e. case of empty DGEresult)
    #stop with an error if item is mismatched with DGEresult dim
    result <- FALSE
    switch(.type[[itemType]],
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
           }

    ) #end switch
    result <- TRUE
    return(result)
}

### addItem
addItem <- function(x, ...) UseMethod("addItem")
addItem.default <- function(dgeResult, ...){
    warning(paste("addItem does not know how to handle object of class ", 
                  class(dgeResult), 
                  "and can only be used on class DGEresult"))
}
addItem.DGEresult <- function(dgeResult, item, itemName, itemType,
                              overwrite=FALSE){
    #enforce itemType
    if (!itemType %in% names(.type))
        stop(paste("itemType must be one of: ", 
                   paste(names(.type), collapse=", "), sep=""))
    
    #add item and assign it's type
    #refuse to add if it exists already unless overwrite = T
    if (overwrite==FALSE & with(dgeResult$data, exists(itemName)))
        stop('itemName already exists in DGEresult!')
    #confirm dimensions consistent before adding
    else if (.dimensionMatch(dgeResult, item, itemType) == TRUE){
        dgeResult$data[[itemName]] <- item
        dgeResult$type[[itemName]] <- itemType
    }
    
    #add workflow info
    if (!with(dgeResult, exists("workflow")))
        dgeResult[["workflow"]] <- list()

    dgeResult$workflow[["ItemName"]] <- list(
            itemName = itemName,
            type = itemType,
            date = date(),
            funArgs = match.call(),
            version = packageVersion("DGE.Tools")
        )
    return(dgeResult)  
} #addItem


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
    dgeResult$data[itemName] <- NULL
    dgeResult$type[itemName] <- NULL
    if (with (dgeResult, exists("workflow")))
        dgeResult$workflow[[itemName]] <- NULL
    return(dgeResult)
}

### getItem
getItem <- function(x, ...) UseMethod("getItem")
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
dim.default <- function(dgeResult, ...) {
    warning(paste("dim does not know how to handle object of class ",
                  class(dgeResult),
                  "and can only be used on class DGEresult"))
}
dim.DGEresult <- function(dgeResult){
    
    dimension <- c(0,0)
    rowcount <- 0
    colcount <- 0
    rowTypes <- names(.type)[.type %in% c("row", "assay")]
    colTypes <- names(.type)[.type %in% c("col", "assay")]
    
    #look for first row or assay or contrastTop object and get nrow
    idx <- dgeResult[["type"]] %in% rowTypes
    if (sum(idx) > 0) {
        firstItemName <- names(dgeResult$type[idx])[[1]]
        rowcount <- nrow(dgeResult$data[[firstItemName]])
    }
    
    #look for first col or assay element
    idx <- dgeResult[["type"]] %in% colTypes
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

#Return all items of a specified type as a list
getType <- function(x, ...) UseMethod("getType")
getType.default <- function(dgeResult, ...){
    warning(paste("dim does not know how to handle object of class ", 
                  class(dgeResult), 
                  "and can only be used on class DGEresult"))
}
getType.DGEresult <- function(dgeResult, type){
    idx <- dgeResult$type %in% type
    itemList <- dgeResult$data[idx]
}  #dim

### print
print <- function(x, ...) UseMethod("print")
print.default <- function(dgeResult, ...) {
    warning(paste("print does not know how to handle object of class ",
                  class(dgeResult),
                  "and can only be used on class DGEresult"))
}
print.DGEresult <- function(dgeResult, ...)  {
    Dim <- dim(dgeResult)
    itemNames <- paste(names(dgeResult$data), collapse=", ")
    itemTypes <- paste(dgeResult$type, collapse=", ")
    print("Items: ", itemNames, sep="")
    print("ItemTypes: ", itemTypes, sep="")
    print(paste("Dimension: [", Dim[1], ", ", Dim[2], "]", sep=""))

    if (with(dgeResult, exists("workflow"))){
     #collect workflow into a dataframe
     #return the dataframe to let the user decide how to print.
     workflowDF <- do.call(rbind, workflow)
    }
}
