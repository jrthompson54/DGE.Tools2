### DGEresult object definition

### # Constructor function for the class
#Construct and empty DGEresult object
DGEresult <- function(item, itemName, itemType) {
    #initialize an empty DGEresult and optionally add the first item
    object <- list(.dim=c(0,0))
    class(object) <- "DGEresult"
    if (!missing(item) & !missing(itemType))
        object <- addItem.DGEresult(object, item, itemName, itemType)
    return (object)
}

.checkDimension <- function(DGEresult, item, itemType){
    #check the dimensions of item against the DGEresult object dim 
    #set the DGEresult dimension if not set (i.e. case of empty DGEresult)
    #stop with an error if item is mismatched with DGEresult dim
    result <- FALSE
    switch(itemType,
           "row" = {
               if (DGEresult$.dim[1] == 0) DGEresult$.dim[1] <- nrow(item)
               else if (DGEresult$.dim[1] != nrow(item))
                   stop('New row object does not match row dimension of DGEresult object')
           },
           
           "col" = {
               if (DGEresult$.dim[2] == 0) DGEresult$.dim[2] <- nrow(item)
               else if (DGEresult$.dim[2] != nrow(item))
                   stop('Rows in new col object does not match col dimension of DGEresult object')
           },
           
           "assay" = {
               if (DGEresult$.dim[1] == 0) DGEresult$.dim[1] <- nrow(item)
               else if (DGEresult$.dim[1] != nrow(item))
                   stop('New assay object does not match row dimension of DGEresult object')
               if (DGEresult$.dim[2] == 0) DGEresult$.dim[2] <- ncol(item)
               else if (DGEresult$.dim[2] != ncol(item))
                   stop('Rows in new col object does not match col dimension of DGEresult object')
           },
           
           "contrastTop" = {
               if (DGEresult$.dim[1] == 0) DGEresult$.dim[1] <- nrow(item)
               else if (DGEresult$.dim[1] != nrow(item))
                   stop('New row object does not match row dimension of DGEresult object')
           }
           
    ) #end switch
    result <- TRUE
    return(result)
}

addItem.DGEresult <- function(DGEresult, item, itemName, itemType){
    #add item and assign it's type
    #refuse to add if it exists already
    if (with(DGEresult, exists(itemName)))
        stop(paste(itemName, ' already exists in ', DGEresult, '!', sep=""))
    else if (.checkDimension(DGEresult, item, itemType) == TRUE)
       DGEresult[[itemName]] <- list(data=item, type=itemType)
    return(DGEresult)  
}

rmItem.DGEresult <- function(DGEresult, itemName)
    #remove the named item
    DGEresult[[itemName]] <- NULL

getItem.DGEresult <- function(DGEresult, itemName){
    return(DGEresult[ItemName]$data)
}

dim.DGEresult <- function(DGEresult)
    return(DGEresult$.dim)

