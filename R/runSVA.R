### Function runSVA ###
#' Function  runSVA
#'
#' Takes an DGEobj from runVoom and tests for surrogate variables.  Adds a new
#' design matrix to the DGEobj with the surrogate variable columns appended (cbind).
#' runVoom should then be run again with the new design matrix to complete the
#' analysis.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords SVA SLOA
#'
#' @param dgeObj A DGEobj with normalized counts and a DesignMatrix.
#'
#' @return dgeObj The DGEobj is returned containing a new design matrix.
#'
#' @examples
#' MyDgeObj = runSVA (MyDgeObj)
#'
#' @importFrom sva sva
#' @import magrittr
#' @importFrom assertthat assert_that
#' @importFrom DGEobj getItem addItem
#' @importFrom stats model.matrix
#'
#' @export
runSVA<- function(dgeObj, designMatrixName){

  assertthat::assert_that(!missing(dgeObj),
              !missing(designMatrixName),
              class(dgeObj)[[1]] == "DGEobj",
              class(designMatrixName)[[1]] == "character",
              with(dgeObj, exists("design")),
              with(dgeObj, exists(designMatrixName)),
              with(dgeObj, exists(paste(designMatrixName, "_Elist", sep="")))
              )

  #Set up a NullFormula and DesignMatrix
  NullFormula = "~ 1"
  Design = DGEobj::getItem(dgeObj, "design")
  NullDesignMatrix = stats::model.matrix(as.formula(NullFormula), Design)

  log2cpm <- DGEobj::getItem(dgeObj, paste(designMatrixName, "_Elist", sep=""))$E
  designMatrix <- DGEobj::getItem(dgeObj, designMatrixName)
  n.sv <- num.sv(log2cpm, designMatrix, method="leek")
  svobj <- sva::sva(log2cpm, designMatrix, NullDesignMatrix, n.sv=n.sv)
  #pull out the surrogate variables
  sv <- svobj$sv

  if (svobj$n.sv > 0) {

      #give them a colname
      colnames(sv) <- paste("sva", 1:ncol(sv), sep="")

      # Add the SVA colums to the DesignMatrix
      designMatrix_SVA <- cbind (designMatrix, sv)

      #capture output
      FunArgs <- match.call()
      #save the svobj
      saveRDS(svobj, "svobj.RDS")
      dgeObj <- addItem(dgeObj, svobj, paste(designMatrixName, "_svobj", sep=""),
                        "svobj", funArgs=FunArgs,
                        parent=designMatrixName)

      #save the new designMatrix
      dgeObj <- addItem(dgeObj, designMatrix_SVA, paste(designMatrixName, "_sva", sep=""),
                        "designMatrix", funArgs=FunArgs,
                        parent=designMatrixName)

  } else tsmsg ("No Surrogate Variables Found. DGEobj is unchanged.")

  return(dgeObj)
}
