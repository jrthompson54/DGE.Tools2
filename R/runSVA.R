### Function runSVA ###
#' Function  runSVA
#'
#' Takes an DGEobj from runVoom and tests for surrogate variables.  Adds a new
#' design matrix to the DGEobj with the surrogate variable columns appended (cbind).
#' runVoom should then be run again with the new design matrix to complete the
#' analysis.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords SVA surrogate variable DGEobj
#'
#' @param dgeObj A DGEobj with normalized counts and a DesignMatrix.
#' @param designMatrixName The itemName of the design matrix in DGEobj
#' @param method method passed to num.sv. ["leek" or "be"; default = "leek"]
#'
#' @return dgeObj The DGEobj is returned containing a new design matrix and an updated design table.
#'
#' @examples
#' MyDgeObj = runSVA (MyDgeObj)
#'
#' @importFrom sva sva num.sv
#' @import magrittr
#' @importFrom assertthat assert_that
#' @importFrom DGEobj getItem addItem
#' @importFrom stats model.matrix
#'
#' @export
runSVA<- function(dgeObj, designMatrixName, method="leek"){

  assertthat::assert_that(!missing(dgeObj),
              !missing(designMatrixName),
              "DGEobj" %in% class(dgeObj),
              "character" %in% class(designMatrixName),
              tolower(method) %in% c("leek", "be"),
              with(dgeObj, exists("design")),
              with(dgeObj, exists(designMatrixName))
              )

  method <- tolower(method)

  #Set up a NullFormula and DesignMatrix
  NullFormula = "~ 1"
  Design = DGEobj::getItem(dgeObj, "design")
  NullDesignMatrix = stats::model.matrix(as.formula(NullFormula), Design)

  # log2cpm <- DGEobj::getItem(dgeObj, paste(designMatrixName, "_Elist", sep=""))$E
  log2cpm <- convertCounts(dgeObj$counts, unit="cpm", log=TRUE, normalize="tmm")
  designMatrix <- DGEobj::getItem(dgeObj, designMatrixName)
  n.sv <- sva::num.sv(log2cpm, designMatrix, method=method)
  svobj <- sva::sva(log2cpm, designMatrix, NullDesignMatrix, n.sv=n.sv)

  #pull out the surrogate variables
  sv <- svobj$sv

  if (svobj$n.sv > 0) {

      #give them a colname
      colnames(sv) <- paste("sv", 1:ncol(sv), sep="")

      # Add the SVA colums to the DesignMatrix
      designMatrix_SVA <- cbind (designMatrix, sv)

      #capture the function call
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
      #add the SV columns to the Design table
      # svacols <- c(ncol(dgeObj$DiseaseTreatment_sva)-svobj$n.sv+1, ncol(dgeObj_sva$DiseaseTreatment_sva))
      # dgeObj$design <- cbind(dgeObj$design, dgeObj_sva$DiseaseTreatment_sva[,svacols[1]:svacols[2]])
      dgeObj$design <- cbind(dgeObj$design, sv)

  } else tsmsg ("No Surrogate Variables Found. DGEobj is unchanged.")

  return(dgeObj)
}
