### Function plotBCV ###
#' Function  plotBCV
#'
#' Creates an edgeR dispersion plot for RNA-Seq QC purposes.  Takes a counts matrix
#' or DGEList for input.  Dispersion is plotted against AveLogCPM.  
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, dispersion, plot, QC
#'
#' @param x A Counts matrix or DGEList (required)
#' @param designMatrix  A design matrix created limma by model.matrix 
#' @param symbolSize (Default=1)
#' @param symbolShape see
#'   http://www.cookbook-r.com/Graphs/Shapes_and_line_types/ (Default = 1)
#' @param symbolColor Default = "darkblue"
#' @param symbolFill  Default = "darkblue"
#' @param symbolAlpha Transparency for the points. Value from 0 to 1. Smaller
#'   indicate more transparency (Default = 0.3)
#' @param linefitSize (Default = 1)
#' @param linefitColor (Default = "yellow")
#' @param lineFit (Default = NULL) Any type supported by geom_smooth (loess is
#'   recommended). Defaults to NULL because enabling this
#'   overloads my poor 32bit PC.
#' @param rugColor (Default = NULL)  Set to Null to disable rug plots.
#'   Note: printing plots with rug plots can be slow (several minutes on my
#'   meager 32bit pc). So I left this disabled by default.  Assign a valid
#'   color to enable
#' @param rugAlpha Transparency for the rug plots (Default = 0.02)
#' @param ... Extra parameters to pass to edgeR::estimateDisp
#'
#' @return a ggplot object
#'
#' @examples
#'    MyGgplot <- plotBCV (MyDGElist)
#'    MyGgplot <- plotBCV (MyDGEobj)
#'
#' @import assertthat DGEobj
#'
#' @export
plotBCV <- function (x,
                      designMatrix,
                            symbolSize = 1,
                            symbolShape = 1,
                            symbolColor = "darkblue",
                            symbolFill = "darkblue",
                            symbolAlpha = 0.3,
                            linefitSize = 1,
                            linefitColor = "yellow",
                            lineFit = NULL,
                            lineType = 1,
                            lineAlpha = 1,
                            rugColor = NULL,
                            rugAlpha = 0.02,
                            ...) {
# browser()
  #dispersion: expect 0.4  for human, 0.1 for inbreed mouse, near 0 for technical reps.
  #  source("~/R/lib/SubsettableListOfArrays.R")

  assert_that(!missing(x),
              !missing(designMatrix),
              class(designMatrix)[[1]] == "matrix")
    
  if (class(x)[[1]] == "DGEList")
    dgelist <- x %>%
        calcNormFactors %>%
        estimateDisp (design=designMatrix, robust=TRUE, ...)
  else dgelist <- x %>%  #process a counts matrix
      as.matrix %>%
      DGEList %>%
      calcNormFactors %>%
      estimateDisp (design=designMatrix, robust=TRUE, ...) 
  
  plotdata <- data.frame(AveLogCPM=dgelist$AveLogCPM, BCV=sqrt(dgelist$tagwise.dispersion))

  MyBCVPlot <- ggplot(plotdata, aes(x=AveLogCPM, y=BCV)) +
    geom_point (size=symbolSize, shape=symbolShape, fill=symbolFill,
                color=symbolColor, alpha=symbolAlpha)

  if (!is.null(lineFit)){
      MyBCVPlot <- MyBCVPlot +
      geom_smooth(method=lineFit, size=linefitSize, color=linefitColor,
                linetype=lineType, alpha=lineAlpha)
  }

  if (!is.null(rugColor)){
      MyBCVPlot = MyBCVPlot + geom_rug(data=plotdata, inherit.aes=FALSE,
                                         color = rugColor,
                                         alpha=rugAlpha,
                                         show.legend=FALSE,
                                         aes(x=AveLogCPM, y=BCV))
  }

  MyBCVPlot <- MyBCVPlot +
    ggtitle ("EdgeR BCV Plot") +
    theme_grey()

  return(MyBCVPlot)
}
