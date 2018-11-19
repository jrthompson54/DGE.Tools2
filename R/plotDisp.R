### Function plotDisp ###
#' Function  plotDisp
#'
#' Creates an edgeR dispersion plot for RNA-Seq QC purposes.  Takes a counts matrix
#' or DGEList for input.  Dispersion is plotted against AveLogCPM.  Optionally,
#' you can plot Biological Coefficient of Variation instead (BCV is the sqrt of
#' dispersion).
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, dispersion, plot, QC
#'
#' @param x A Counts matrix or DGEList (required)
#' @param designMatrix  A design matrix created stats::model.matrix
#' @param plotType One of "dispersion" or "BCV" (default = "dispersion")
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
#'    MyGgplot <- plotDisp (MyDGElist)
#'    MyGgplot <- plotDisp (MyDGEobj)
#'
#' @importFrom assertthat assert_that
#' @importFrom edgeR calcNormFactors estimateDisp DGEList
#'
#' @export
plotDisp <- function (x,
                      designMatrix,
                      plotType = "dispersion",
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
  #dispersion: expect 0.4 coef of var for human, 0.1 for inbreed mouse, near 0 for technical reps.
  #  source("~/R/lib/SubsettableListOfArrays.R")

  assertthat::assert_that(!missing(x),
              !missing(designMatrix),
              class(designMatrix)[[1]] == "matrix")

  if (class(x)[[1]] == "DGEList")
    dgelist <- x %>%
        edgeR::calcNormFactors() %>%
        edgeR::estimateDisp (design=designMatrix, robust=TRUE, ...)
  else dgelist <- x %>%  #process a counts matrix
      as.matrix %>%
      edgeR::DGEList() %>%
      edgeR::calcNormFactors() %>%
      edgeR::estimateDisp (design=designMatrix, robust=TRUE, ...)

  if (tolower(plotType == "dispersion"))
      plotdata <- data.frame(AveLogCPM=dgelist$AveLogCPM, Dispersion=dgelist$tagwise.dispersion)
  else
      plotdata <- data.frame(AveLogCPM=dgelist$AveLogCPM, Dispersion=sqrt(dgelist$tagwise.dispersion))

  MyDispPlot <- ggplot(plotdata, aes(x=AveLogCPM, y=Dispersion)) +
      geom_point (size=symbolSize, shape=symbolShape, fill=symbolFill,
                  color=symbolColor, alpha=symbolAlpha)

  if (!is.null(lineFit)){
   MyDispPlot <- MyDispPlot +
      geom_smooth(method=lineFit, size=linefitSize, color=linefitColor,
                linetype=lineType, alpha=lineAlpha)
  }

  if (!is.null(rugColor)){
    MyDispPlot = MyDispPlot + geom_rug(data=plotdata, inherit.aes=FALSE,
                                         color = rugColor,
                                         alpha=rugAlpha,
                                         show.legend=FALSE,
                                         aes(x=AveLogCPM, y=Dispersion))
  }

  MyDispPlot <- MyDispPlot +
    ggtitle ("EdgeR Dispersion Plot") +
    expand_limits(x=0, y=0) +
    theme_grey()

  if (tolower(plotType) == "bcv")
      MyDispPlot <- MyDispPlot +
                    ylab("BCV") +
                    ggtitle("EdgeR BCV Plot")

  return(MyDispPlot)
}
