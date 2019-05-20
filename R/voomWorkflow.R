### Function voomWorkflow ###
#' Function voomWorkflow
#'
#' This function runs several steps in the RNA-Seq pipeline in one fell swoop.
#' It supports duplicateCorrelation if you provide a blocking vector.  It
#' includes low intensity filtering, filtering for protein coding genes and
#' filters out zero effective length genes (genes shorter than library size).
#' Then it runs TMM normalization, voomWithQualityWeights, lmFit and eBayes.
#'
#' To incorporate SVA analysis, you need to run SVA first and add the SVA
#' variables to your design table. Then you can proceed with this function.
#'
#' After running this step. You define your contrasts and execute runContrast to
#' complete DGE calculations.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, DGEobj, limma voom
#'
#' @param dgeObj  A class dgeObj with counts, gene annotation and sample
#'   annotation
#' @param formula A text representation of the formula you want to use (clas
#'   character not formula)
#' @param projectName This should be the project name from Xpress or Omicsoft
#' @param designMatrixName User defined name for the design matrix
#' @param dupcorBlock A blocking vector to define which samples belong to the
#'   same subject to be used with the duplicateCorrelation function.
#' @param outputPath Where to send output plots
#' @param annotationFile Text file of key=value pairs to populate DGEobj
#'   attributes (optional but highly advised)
#' @param proteinCodingOnly Set to TRUE to keep only protein coding genes
#'   (default = TRUE)
#' @param sampleFraction Fraction of samples that must meet intensity thresholds
#'   to keep a gene (Default = 0.5)
#' @param ... Additional named arguments passed to the low intensity filter
#'   function to define the desired intensity filter type (see ?lowIntFilter). Settable
#'   arguments for low intensity filtering are: fracThreshold (default = 0.5),
#'   countThreshold, fpkThreshold, zfpkmThreshold, tpmThreshold.  You can use
#'   countThreshold plus one other argument.  If no arguments supplied here, the
#'   following defaults apply (fracThreshold=0.5, countThreshold=10,
#'   fpkThreshold=5).  To disable intensity filtering use: sampleFraction = 0.
#'
#' @return A DGEobj with analysis results added
#'
#' @examples
#'
#'    MyDgeObj <- voomWorkflow(MyDgeObj,
#'                             formula = "~ 0 + ReplicateGroup",
#'                             projectName = "MyProjectName",
#'                             designMatrixName = "ReplicateGroup",
#'                             annotationFile = "MyProjectName.txt",
#'                             proteinCodingOnly = TRUE,
#'                             )
#'
#' @import magrittr DGEobj
#' @importFrom assertthat assert_that
#'
#' @export
voomWorkflow <- function(dgeObj,
                         formula,
                         projectName,
                         designMatrixName,
                         dupcorBlock,
                         outputPath = "./",
                         annotationFile,
                         proteinCodingOnly = FALSE,
                         sampleFraction=0.5,
                         ...){

  assertthat::assert_that(!missing(dgeObj),
                          !missing(formula),
                          !missing(projectName),
                          !missing(designMatrixName),
                          "DGEobj" %in% class(dgeObj),
                          "character" %in% class(formula),
                          "character" %in% class(projectName),
                          "character" %in% class(designMatrixName))

  ellipsisArgs <- list(...)

  RDSname <- file.path(outputPath, str_c(projectName, ".RDS"))

  #add project metadata
  if (!missing(annotationFile))
    dgeObj <-  annotateDGEobj(dgeObj, regfile=annotationFile)

  #check for and create output folder if needed
  if (outputPath != "./" & !file.exists(outputPath))
    dir.create(file.path(outputPath))

  #Filter out genes with any zero efflength (only needed for data from Xpress)
  suppressMessages(el <- getItem(dgeObj, "effectiveLength"))
  if (!is.null(el)){  #only Xpress data contains an effectiveLength item
    rowmin <- apply(el, 1, min)
    idx <- rowmin > 0
    dgeObj <- dgeObj[idx,]
  }


  #trap for too many intensity filter arguments
  if (length(ellipsisArgs) > 2) stop("Only 1-2 intensity filtering args allowed")
  if (length(ellipsisArgs) == 2 & !"countThreshold" %in% names(ellipsisArgs))
    stop("Two intensity filtering args were supplied, one of them must be countThreshold")
  #now we should have 1-2 filtering args

  #intensity filtering  args from ellipsis
  if ("fracThreshold" %in% names(ellipsisArgs))
    fracThreshold <- ellipsisArgs$fracThreshold
  if ("fpkThreshold" %in% names(ellipsisArgs))
    fpkThreshold <- ellipsisArgs$fpkThreshold
  if ("countThreshold" %in% names(ellipsisArgs))
    countThreshold <- ellipsisArgs$countThreshold
  if ("zfpkmThreshold" %in% names(ellipsisArgs))
    zfpkmThreshold <- ellipsisArgs$zfpkmThreshold
  if ("tpmThreshold" %in% names(ellipsisArgs))
    tpmThreshold <- ellipsisArgs$tpmThreshold

  #construct filtering command
  cmd <- ("dgeObj <- lowIntFilter(dgeObj, sampleFraction=sampleFraction, ")
  for (i in 1:length(ellipsisArgs)){
    cmd <- str_c(cmd, names(ellipsisArgs)[i], " = ", ellipsisArgs[[i]])
    if (i < length(ellipsisArgs)) cmd <- str_c(cmd, ", ")
  }
  cmd <- str_c(cmd, ")")
  eval(parse(text=cmd))

  #keep only protein coding genes
  if (proteinCodingOnly == TRUE){
    idx <- NULL
    if ("Source" %in% colnames(dgeObj$geneData)){ #Omicsoft field is Source
      idx <- dgeObj$geneData$Source == "protein_coding"
    } else if ("Source" %in% colnames(dgeObj$isoformData)){
      idx <- dgeObj$isoformData$Source == "protein_coding"
    }

    if (is.null(idx) == FALSE){
      #Xpress sometines has NA in the Source/Biotype field; need to convert those to FALSE
      idx [is.na(idx)] <- FALSE
      dgeObj <- dgeObj[idx,]
    }
  }

  #set up and save the design matrix
  designMatrix <- model.matrix (as.formula(formula), getItem(dgeObj, "design"))
  #clean up problem characters in colnames
  colnames(designMatrix) <- make.names(colnames(designMatrix))
  #capture the formula as an attribute
  designMatrix <- setAttributes(designMatrix, list(formula=formula))
  #save the modified designMatrix to the DGEobj
  dgeObj <- addItem(dgeObj, item=designMatrix,
                    itemName=designMatrixName,
                    itemType="designMatrix",
                    parent="design")


  ### run DGE calculations
  #normalize
  dgeObj <- runEdgeRNorm(dgeObj, plotFile=file.path(outputPath, str_c("TMM_Norm.Factors.PNG")), normMethod = "TMM")

  message("Running Voom... this takes some time...")

  #voom/lmfit
  if (missing(dupcorBlock) || is.null(dupcorBlock)){
    dgeObj <- runVoom(dgeObj, designMatrixName,
                      qualityWeights = TRUE,
                      mvPlot = TRUE)
  } else {
    dgeObj <- runVoom(dgeObj, designMatrixName,
                      qualityWeights = TRUE,
                      dupcorBlock = dupcorBlock,
                      mvPlot = TRUE)
  }

}


