### Function voomWorkflow ###
#' Function voomWorkflow
#'
#' This function runs several steps in the RNA-Seq pipeline in one fell swoop.
#' It supports duplicateCorrelation if you provide a blocking vector.  It does
#' not directly support SVA.
#'
#' For SVA, you need to run SVA first and add the SVA variables to your design
#' table. Then you can proceed with this function.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, DGEobj
#'
#' @param dgeObj  A class dgeObj with counts, gene annotation and sample annotation
#' @param formula A text representation of the formula you want to use
#' @param projectName This should be the project name from Xpress or Omicsoft
#' @param designMatrixName User defined name for the design matrix
#' @param fracThreshold Fraction of samples that must meet intensity thresholds to keep a gene (Default = 0.5)
#' @param outputPath Where to send output plots
#' @param annotationFile File of key=value pairs to populate required attributes (optional but highly advised)
#'
#' @return A DGEobj with a modified object definition embedded.
#'
#' @examples
#'
#'    MyDgeObj <- newType(MyDgeObj,
#'
#' @import magrittr DGEobj
#' @importFrom assertthat assert_that
#'
#' @export
voomWorkflow <- function(dgeObj,
                         formula,
                         projectName,
                         designMatrixName,
                         fracThreshold = 0.5,
                         outputPath = "./",
                         annotationFile){

  assertthat::assert_that(!missing(dgeObj),
                          !missing(formula),
                          !missing(projectName),
                          !missing(designMatrixName),
                          "DGEobj" %in% class(dgeObj),
                          "character" %in% class(formula),
                          "character" %in% class(projectName),
                          "character" %in% class(designMatrixName))

  RDSname <- file.path(outputPath, str_c(projectName, ".RDS"))

  #add project metadata
  if (!missing(annotationFile))
    dgeObj <-  annotateDGEobj(dgeObj, regfile=annotationFile)

  #check for and create output folder if needed
  if (outputPath != "./" & !file.exists(outputPath))
    dir.create(file.path(outputPath))

  #Filter out genes with any zero efflength (only needed for data from Xpress)
  el <- getItem(dgeObj, "effectiveLength")
  if (!is.null(el)){  #only Xpress data contains an effectiveLength item
    rowmin <- apply(el, 1, min)
    idx <- rowmin > 0
    dgeObj <- dgeObj[idx,]
  }

  #Low Intensity Filter
  dgeObj <- lowIntFilter(dgeObj, zfpkmThreshold = -3, countThreshold = 10, sampleFraction=fracThreshold)

  #set up and save the design matrix
  designMatrix <- model.matrix (as.formula(formula), getItem(dgeObj, "design"))
  #clean up problem characters in colnames
  colnames(designMatrix) <- make.names(colnames(designMatrix))
  #capture the formula as an attribute
  designMatrix <- setAttributes(designMatrix, list(formula=formula))
  #save the modified designMatrix
  dgeObj <- addItem(dgeObj, item=designMatrix,
                    itemName=designMatrixName,
                    itemType="designMatrix",
                    parent="design", overwrite=TRUE)


  ### run DGE calculations
  #normalize
  dgeObj <- runEdgeRNorm(dgeObj, plotFile=file.path(outputPath, str_c("TMM_Norm.Factors.PNG")), normMethod = "TMM")  #TMM is the default

  message("Running Voom")

  #voom/lmfit
  if (is.null(dupcorBlock)){
    dgeObj <- runVoom(dgeObj, designMatrixName,
                      qualityWeights = TRUE,
                      mvPlot = FALSE)
  } else {
    dgeObj <- runVoom(dgeObj, designMatrixName,
                      qualityWeights = TRUE,
                      dupcorBlock = dupcorBlock)
  }

}


