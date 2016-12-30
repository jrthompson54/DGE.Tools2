### Function runEdgeRNorm ###
#' Function  runEdgeRNorm (edgeR normlization )
#'
#' Returns a SubsettableListOfArrays (SLOA) containing the Fit and associated
#' output from the standard RNA-Seq using edgeR TMM normalization
#' followed by limma voom and lmfit.  The SLOA contains the Fit object
#' and various associated data  use "print(SLOA)" to examine the contents.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID
#'
#' @param dat A DGEobj or RSE object data structure containing counts and design
#'    data
#' @param normMethod One of "TMM", "RLE", "upperquartile" or "none". Default = "TMM"
#' @param plotFile Enable a Barplot of the norm.factors produced (Default = "Norm.Factors.PNG")
#'   Set to NULL to disable the plot.
#' @param plotLabels Sample labels for the plot. Length must equal the number of
#'   samples. (Default = NULL;  sample number will be displayed)
#'
#' @return A DGEObj with a normalized DGEList added.
#'
#' @examples
#' MyDgeObj <- runEdgeRNorm (RSE)
#' 
#' MyDgeObj <- runEdgeRNorm (MyDgeObj)
#'
#' @import edgeR magrittr limma SummarizedExperiment DGEobj assertthat
#'
#' @export
runEdgeRNorm <- function(dat, NormMethod="TMM", 
                         plotFile="TMM_Norm.Factors.PNG",
                         plotLabels = NULL){

    funArgs <- match.call()
    
    datClass <- class(dat)[[1]]
    assert_that ((datClass == "RangedSummarizedExperiment") |
                (datClass == "DGEobj"))
   
#create the DGEobj if needed  
  if (datClass == "RangedSummarizedExperiment")
      dgeObj <- as.DGEO(dat)
  else dgeObj <- dat

  #convert counts to a matrix
  CountsMatrix = as.matrix(getType(dgeObj, "counts"))

  #Now ready to normalize counts
  MyDGElist = CountsMatrix %>%
    DGEList %>%  #edgeR
    calcNormFactors (method = NormMethod)   
  
  #capture the DGEList
  custAttr <- list(normalization="TMM")
  dgeObj <- addItem(dgeObj, MyDGElist, "DGEList",  "DGEList", 
                    funArgs=funArgs,
                    custAttr=custAttr)

  #plot the Norm factors
  if (!is.null(plotLabels)  && length(plotLabels == ncol(dgeObj))){
    x = plotLabels
    angle = 45
  } else {
    x = 1:ncol(dgeObj)
    angle = 0
  }

  if (!is.null(plotFile)){ #barplot of norm.factors
    df <- data.frame(x = factor(x),
                     Norm.Factors = MyDGElist$samples$norm.factors)
    nfplot <- ggplot (df, aes(x=x, y=Norm.Factors)) +
                geom_bar (stat="identity", color = "dodgerblue4",
                          fill = "dodgerblue3", width = 0.7) +
                geom_hline(yintercept = 1.0, color = "red") +
                xlab("Samples") +
                ylab("Norm Factors") +
                ggtitle ("Normalization Factors") +
                theme_bw(12) +
                theme(axis.text.x = element_text(angle=angle, hjust = 1.0))

    print(nfplot)
    printAndSave (nfplot, filename = plotFile, printPlot=FALSE)

  }

  return(dgeObj)
}
