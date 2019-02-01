### Function QCplots ###
#' Function QCplots
#'
#' This function takes a dataframe with metric names in col1 and samples in col2-n.  Each row is a different metric.
#' It then generates the qc pltos as defined by other arguments.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords barplot lineplot gplot2 logratio confidence intervals contrasts
#'
#' @param qcdata A dataframe or tibble with metric names in the first column and
#'   samples in columns 2-n.  Each row is a different QC metric. This matches
#'   the Omicsoft RNA-Seq.QCMetrics.Table.txt format.
#' @param metricNames A list of metrics to plot.  Values must exist in column 1
#'   of the data frame.
#' @param plotType One of "bar", "point", "pointline".  If you want a different
#'   plottype for each metric, pass a list of plotTypes with length equal to
#'   length(metricNames)
#' @param barColor Color for the bar outline (default = "dodgerblue4")
#' @param barFill Color for the bar area (default = "dodgerblue3")
#' @param barSize set the bar size (thickness of each bar perimeter; default =
#'   0.1)
#' @param barWidth set the bar width (Default = 0.8)
#' @param barAlpha Transparency for the bar layer (Default = 1)
#' @param pointColor Color for the point layer (Default = "grey30")
#' @param pointFill Fill color for the point layer (Default = "dodgerblue4")
#' @param pointShape Shape for the point layer (Default = 21; fillable circle)
#' @param pointAlpha Transparency for the box layer (Default = 1)
#' @param pointSize Size of the points (Default = 4)
#' @param lineColor Color of the line fit (Default = "dodgerblue4")
#' @param lineSize Size of the line fit (Default = 1)
#' @param lineFit Type of fit to use.  One of c("auto", "lm", "glm", "gam",
#'   "loess"). (Default = "loess")
#' @param lineType One of c("solid", "dashed", "dotted", "dotdash", "longdash",
#'   "twodash"). (Default = "solid")
#' @param xAngle Angle to set the sample labels on the Xaxis (Default =  45;
#'   Range = 0-90)
#' @param debug Turn on debug mode (default = FALSE)
#'
#' @return ggplot object if one plot is specified.  A list of ggplot objects if 2 or more metrics specified
#'
#' @examples
#'
#'
#' @import ggplot2 magrittr
#' @importFrom assertthat assert_that
#' @importFrom stringr str_c
#'
#' @export
QCplots <- function(qcdata,
                    metricNames,
                    plotType = "bar",
                      barColor = "dodgerblue4",
                      barFill = "dodgerblue3",
                      barSize = 0.1,
                      barAlpha = 1,
                      barWidth = 0.9,
                      pointColor = "grey30",
                      pointFill = "dodgerblue4",
                      pointShape = 21, #fillable circle
                      pointAlpha = 1,
                      pointSize = 2,
                      lineColor = "dodgerblue4",
                      lineSize = 1,
                      lineType = "solid",
                      lineFit = "loess",
                      xAngle = 45,
                      debug = FALSE
                      )
{
  assertthat::assert_that("data.frame" %in% class(qcdata),
                          tolower(plotType) %in% c("bar", "point", "pointline"),
                          all(metricNames %in% qcdata[,1])
  )

  #convert first col to rownames
  qcdata %<>% column_to_rownames(var=colnames(qcdata)[1]) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var="Sample")

  #if only one plotType, apply to all plots
  if(length(plotType) ==1) plotType <- rep(plotType, length(metricNames))

  plots <- list()
  for (metric in metricNames){
    idx <- metric == metricNames
    plot_type <- plotType[idx]

    p <- ggplot(qcdata, aes(x=Sample, y=!!metric))

    p <- switch(tolower(plot_type),
           bar = p + geom_bar(color=barColor, fill=barFill, alpha=barAlpha, width=barWidth),
           point = p + geom_point(color=pointColor, fill=pointFill, shape=pointShape, alpha=pointAlpha),
           pointline = {p +
               geom_point(color=pointColor, fill=pointFill, shape=pointShape, alpha=pointAlpha) +
               geom_smooth(method="lm", color=lineColor, size=lineSize, linetype=lineType, fit=lineFit)}
    )
    plots <- list(plots, p)
  }
  if (length(plots == 1)){
    return(plots[[1]])
  } else {
    return(plots)
  }

}

