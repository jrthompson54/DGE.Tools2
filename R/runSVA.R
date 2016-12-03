### Function runSVA ###
#' Function  runSVA
#'
#' Takes an SubsettableListOfArrays (SLOA) object from runVoom and returns
#' the SLOA object after running SVA and adding the surrogate variables it
#' discovers to the DesignMatrix in the SLOA.  The resulting SLOA object is
#' ready to be reprocessed by runVoom again to complete the analysis.
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
#' @import sva SummarizedExperiment magrittr
#'
#' @export
runSVA<- function(MySLOA){

  #Version Info
  FVersion = "runSVA : 28Dec2015"

  if (!exists("SLOA")) {
    stop("You need to source SubsettableListOfArrays.R before executing runSVA")
  }

  #Set up a NullFormula and DesignMatrix
  NullFormula = "~ 1"
  Design = MySLOA$colData
  NullDesignMatrix = model.matrix(as.formula(NullFormula), Design)

  n.sv <- num.sv(MySLOA$Elist$E, MySLOA$DesignMatrix, method="leek")
  svobj <- sva(MySLOA$Elist$E, MySLOA$DesignMatrix, NullDesignMatrix, n.sv=n.sv)

  # Add the SVA columns to the Design table
  NewDesign <- MySLOA$colData %>% as.data.frame
  NewDesign %<>% cbind(svobj$sv)
  MySLOA$colData <- DataFrame(NewDesign)

  # Add the SVA colums to the DesignMatrix
  MySLOA$DesignMatrix %<>% cbind (svobj$sv)

  #add the svobj object to the MySLOA
  MySLOA$svobj <- svobj
  MySLOA$workflowRecord$SVA <- paste (date(), " : v", packageVersion("DGE.Tools"), sep="")
  MySLOA$workflowRecord$Date = date()

  #delete the old Fit and Residuals to avoid confusion
  MySLOA$Fit <- NULL
  MySLOA$Residuals <- NULL

  return(MySLOA)
}
