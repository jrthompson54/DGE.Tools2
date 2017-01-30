### Function runSVA ###
#' Function  runSVA
#'
#' Takes an DGEobj from runVoom and tests for surrogate variables.  Adds a new
#' design matrix to the DGEobj with the surrogate varaible columns appended. 
#' runVoom should then be run again with the new design matrix to complete the
#' analysis.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords SVA SLOA
#'
#' @param MySLOA An SLOA object with normalized counts and a DesignMatrix.
#'
#' @return MySLOA An SLOA object containing data ready for DGE analysis (runVoom)
#'
#' @examples
#' MySLOA = runSVA (MySLOA)
#'
#' @import sva magrittr assertthat
#'
#' @export
runSVA<- function(dgeObj, designMatrixName){

  assert_that(!missing(dgeObj),
              !missing(designMatrixName),
              class(dgeObj)[[1]] == "DGEobj",
              class(designMatrixName)[[1]] == "character",
              with(dgeObj, exists("design")),
              with(dgeObj, exists(designMatrixName)),
              with(dgeObj, exists(paste(designMatrixName, "_Elist", sep="")))
              )

  #Set up a NullFormula and DesignMatrix
  NullFormula = "~ 1"
  Design = getItem(dgeObj, "design")
  NullDesignMatrix = model.matrix(as.formula(NullFormula), Design)
  
  log2cpm <- getItem(dgeObj, paste(designMatrixName, "_Elist", sep=""))$E
  designMatrix <- getItem(dgeObj, designMatrixName)
  n.sv <- num.sv(log2cpm, designMatrix, method="leek")
  svobj <- sva(log2cpm, designMatrix, NullDesignMatrix, n.sv=n.sv)

  # Add the SVA colums to the DesignMatrix
  designMatrix_SVA <- cbind (designMatrix, svobj$sv)
  
  #capture output
  FunArgs <- match.call()
  #save the svobj
  saveRDS(svobj, "svobj.RDS")
  dgeObj <- addItem(dgeObj, svobj, paste(designMatrixName, "_svobj", sep=""),
                    "svobj", funArgs=FunArgs,
                    custAttr=list(parent=designMatrixName))
  #save the new designMatrix
  dgeObj <- addItem(dgeObj, designMatrix_SVA, paste(designMatrixName, "_sva", sep=""),
                    "designMatrix", funArgs=FunArgs,
                    custAttr=list(parent=designMatrixName))

  return(dgeObj)
}
