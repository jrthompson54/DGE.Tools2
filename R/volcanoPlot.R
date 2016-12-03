### Function volcanoPlot ###
#' Deluxe Volcano Plots
#'
#' A volcano plot shows Log Ratio data on the X axis and Negative Log Pvalues (NLP) on the Y axis.
#' This function is intended to show the volcano plot from a dataframe
#' created by topTable or topTreat.  Properly normalized data will generally be
#' centered around LogRatio = 0.
#'
#' By default, the plot places "logFC" on the X axis and Log10 of the "P.Value" on the Y axis.
#' By default, a reference vertical line is drawn at LogRatio = 0 on the X axis.
#' Optionally, additional reference lines will be drawn at +/- a user supplied
#' logratio threshold.
#' The points are color coded using both the significance and fold-change thresholds supplied by the user.
#' By default, the P.Value field is used with a threshold of 0.01 to color code the points and fold-change
#' threshold of +/- 1.5X.
#'
#' The volcanoPlot has also been implemented
#' with a scalable theme that makes it easy to scale font sizes larger
#' or smaller for PPT or Knitr output.  Just add baseFont(n),
#' where n equals the desired base font size (12 works well for knitr, 18 or 24 works well
#' for PPT).  e.g. MyPlot + baseFont(18).
#'
#' \strong{Data Structure for the input dataframe:}
#'
#' The defaults are set for dataframes produced by topTable and topTreat.  The columns named "logFC"
#' and "P.Value" are used by default to accommodate the column
#' names used in topTable/topTreat dataframes.  Any other dataframe
#' can be used with foldchange, intensity and significance measures.  In this case you must use
#' the appropriate arguments to define the column names to use. By default, the
#' column names will be used for the axis labels.  You can overide the default xylabels
#' the xlab and ylab arguments.
#'
#' A significance measure (which defaults to P.Value <= 0.01) and LogRatio
#' threshold are used to color code genes that
#' are significantly increased or decreased.
#' Use the appropriate arguments to use an FDR measure instead of pvalue.
#'
#' Sensible defaults are chosen for symbols (Size, Shape, Color and Fill).  Optionally,
#' there are arguments to allow you to fiddle with these settings.  A length of 3 is
#' required for these arguments which applies the attributes in this order:
#' Increased, NoChange, Decreased.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords compare plot ggplot2 logratio scatterplot
#'
#' @param df A dataframe with LogRatio and LogIntensity columns and optionally a
#'   pvalue or FDR column.
#' @param logRatioCol name of the LogRatio column (Default = "logFC")
#' @param logIntCol name of the LogIntensity column (Default = "AveExpr")
#' @param pvalCol name of the pvalue or FDR column (Default = "P.Value")
#' @param xlab X axis label (default the LogIntensity column name)
#' @param ylab Y axis label (default the LogRatio column name)
#' @param title Plot title (optional)
#' @param pthreshold Used to color points (default = 0.01)
#' @param pthresholdLine Color for a horizontal line at the pthreshold (Default
#'   = NULL (disabled))
#' @param symbolSize Size of symbols for Up, no change and Down. default = c(4,
#'   3.99, 4); Note: you cannot choose the exact same size for all three.  But
#'   you can use decimal values.
#' @param symbolShape Shape of the symbols for Up, no change and Down; Default =
#'   c(21, 1, 21) (1 = open circle, 21 = fillable open circle); Note: you cannot
#'   choose the same symbol shape for all three. See
#'   \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Up, NoChange, Down); default = c("black", "grey25",
#'   "grey0") See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#'   Note: You cannot duplicate colors.
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25
#'   are fillable. This will have no effect on other symbols. Default =
#'   c("red3", "grey25", "deepskyblue4") Note: You cannot duplicate colors.
#' @param scaleByIntensity Set point size as a function of Intensity (Default =
#'   TRUE). Takes logIntCol and sets a floor and ceiling at 0 and 10, then uses
#'   this value to size the points.
#' @param alpha Controls the transparency of the plotted points (range: 0-1;
#'   default = 0.7)
#' @param foldChangeLines Position of reference vertical lines for fold change
#'   (Default = log2(1.5); NULL disables)
#' @param legendPosition One of "top", "bottom", "left", "right", "ne", "se",
#'   "nw", "sw", NULL. top/bottom/left/right place the legend outside the
#'   figure.  ne/se/nw/sw place the figure inside the figure. NULL disables the
#'   legend. Default = "right"
#' @param rugColor Specify color for rug density plot along xy axes. Set to NULL
#'   to disable rug layer. (Default=NULL)
#' @param rugAlpha Sets the transparency for the rug layer.  An alpha <1.0 takes
#'   a long time to draw (~10 min on 32bit Win7).  So I set the default at 1.0,
#'   but for a final plot try 0.1 or 0.2 for alpha for a more informative rug.
#' @param baseFontSize The smallest size font in the figure in points. Default =
#'   12
#' @param themeStyle "bw" or "grey" which correspond to bwTheme or greyTheme
#'   respectively. Default = bw"
#' @param refLineThickness Set the thickness for all reference lines (Default =
#'   1)
#'
#' @return ggplot object
#'
#' @examples
#'    Simple plot with custom title(df is a topTable dataframe)
#'
#'    MyPlot = volcanoPlot(df, title = "Plot Title")
#'
#'    Some options with a custom datafile
#'
#'    MyPlot = volcanoPlot (df, pthreshold = 0.1,
#'      logRatioCol = "Log2ratio",
#'      logIntCol = "AverageIntensity",
#'      pvalCol = "BHFDR",
#'      xlab = "Log2 Ratio", ylab = "Log10 BHFDR",
#'      title = "Profile Plot Title",
#'      referenceLine = "blue",
#'      legendPosition="ne")
#'
#' @import ggplot2 magrittr dplyr
#'
#' @export
volcanoPlot <- function(df,
                        logRatioCol = "logFC",
                        logIntCol = "AveExpr",
                        pvalCol = "P.Value",
                        pthreshold=0.01,
                        xlab=NULL, ylab=NULL, title=NULL,
                        symbolSize = c(4, 3.999, 4),
                        symbolShape = c(21, 1, 21),
                        symbolColor = c("black", "grey25", "grey0"),
                        symbolFill = c("red3", "grey25", "deepskyblue4"),
                        alpha = 0.5,
                        sizeByIntensity = TRUE,
                        pthresholdLine = NULL,
                        foldChangeLines = log2(1.5),
                        refLineThickness = 1,
                        legendPosition = "right",
                        rugColor = NULL,
                        rugAlpha = 1.0,
                        baseFontSize = 12,
                        themeStyle = "grey"
                        ) {

#   basetheme = theme(  #base size is the axis tick mark labels; other elements are scaled
#     axis.text.x = element_text(size=rel(1.0)),
#     axis.text.y = element_text(size=rel(1.0)),
#     axis.title.x = element_text(size=rel(1.25), vjust = 0.5, hjust=0.5, color="black"),
#     axis.title.y = element_text(size=rel(1.25), vjust = 0.5, hjust=0.5, color="black"),
#     plot.title = element_text(face="bold", size = rel(1.5)),
#     legend.text = element_text(colour="Black", size=rel(1.0)),
#     legend.title = element_text(colour="Black", size=rel(1.2)),
#     strip.text.x = element_text(size=rel(0.8)),
#     strip.text.y = element_text(size=rel(0.8))
#   )

  #argument checks
  #
  # Make sure specified columns exist
  if (!logRatioCol %in% colnames(df)) {
    stop("LogRatio Ccolumn not found.")
  }
  if (!logIntCol %in% colnames(df)) {
    stop("LogIntensity column not found.")
  }
  if (!pvalCol %in% colnames(df)) {
    stop("Significance measure column not found.")
  }

  #symbol parameters must all be length=3
  if (!length(symbolSize)==3 || !length(symbolShape)==3 ||
      !length(symbolColor)==3 || !length(symbolFill)==3){
    stop("symbol arguments must be length=3 (symbolSize, SymbolShape, SymbolColor, symbolFill)")
  }

#   #all symbolSize values cannot be the same or ggplot throws an error.
#   if (symbolSize[1] == symbolSize[2]) {
#     symbolSize[2] = symbolSize[2] - 0.01
#   }
#
#   #Check for unique colors (all three colors must be unique for the legend to work properly,
#   #even if the fill color is being applied to a non-fillable symbol
#   if (length(unique(symbolColor)) < 3){
#     stop("All 3 symbol colors must be unique")
#   }
#   if (length(unique(symbolFill)) < 3){
#     stop("All 3 symbol fill colors must be unique")
#   }
#   if (length(unique(symbolShape)) < 2){
#     stop("You need to use at least 2 different symbol shapes")
#   }

  if (sizeByIntensity==TRUE){
    #create a column to support sizeByIntensity
    df$LogInt = df[[logIntCol]]
    #set a floor and a ceiling
    df$LogInt[df$LogInt<0] = 0
    df$LogInt[df$LogInt>10] = 10
  }

  names(symbolShape) = c("Increased", "No Change", "Decreased")
  names(symbolSize) = c("Increased", "No Change", "Decreased")
  names(symbolColor) = c("Increased", "No Change", "Decreased")
  names(symbolFill) = c("Increased", "No Change", "Decreased")

  ssc = data.frame(group = c("Increased", "No Change", "Decreased"),
                   symbolShape = symbolShape,
                   symbolSize = symbolSize,
                   symbolColor = symbolColor,
                   symbolFill = symbolFill,
                   order = c(1,3,2),
                   stringsAsFactors = FALSE) %>% arrange(order)
                  #setting the order defines the legend order

  #columns to plot
  #capture the labels from the colname
  xlabel = logRatioCol
  ylabel = paste ("-log10(", pvalCol, ")", sep="")
  #now make the columnames suitable for use with aes_string
  x = make.names(colnames(df)[colnames(df)==logRatioCol])
  colnames(df)[colnames(df)==logRatioCol] = make.names(colnames(df)[colnames(df)==logRatioCol])
  #make a log10significance column and make that the y column
  df$NegativeLogP = -log10(df[,pvalCol])
  y = "NegativeLogP"

  ### DELUXE PLOT: plot groups in different colors/shapes

  #Let's plot the subsets
  DEup = df[[pvalCol]] <= pthreshold & df[[logRatioCol]] > 0
  DEdn = df[[pvalCol]] <= pthreshold & df[[logRatioCol]] < 0
  DEnot = !DEup & !DEdn
  #create group factor column in df
  df$group=NA
  df$group[DEup] = "Increased"
  df$group[DEdn] = "Decreased"
  df$group[DEnot] = "No Change"
  df %<>% left_join(ssc)
  df$group %<>% factor(levels=c("Increased", "Decreased", "No Change"))

  #set an order field to force plotting of NoChange first
  df$order = NA
  df$order[DEup] = 1
  df$order[DEdn] = 1
  df$order[DEnot] = 0

  volcanoPlot <- ggplot (df, aes_string(x=x, y=y)) +
                    aes(shape=group, size=group,
                      color=group, fill=group,
                      order=order) +
      #scale lines tell it to use the actual values, not treat them as factors
      scale_shape_manual(name="Group", guide="legend", labels=ssc$group,
                           values=ssc$symbolShape) +
      scale_size_manual(name="Group", guide="legend", labels=ssc$group,
                          values=ssc$symbolSize) +
      scale_color_manual(name="Group", guide="legend", labels=ssc$group,
                           values=ssc$symbolColor) +
      scale_fill_manual(name="Group", guide="legend", labels=ssc$group,
                          values=ssc$symbolFill) +
      geom_point(alpha=alpha) +
      #box around the legend
      theme(legend.background = element_rect(fill="gray95", size=.5, linetype="dotted"))
###end DELUXE PLOT

  ### Optional Decorations

  if (!is.null(rugColor)){
    volcanoPlot = volcanoPlot + geom_rug(data=df, inherit.aes=FALSE,
                                         color = rugColor,
                                         alpha=rugAlpha,
                                         show.legend=FALSE,
                                         aes_string(x=x, y=y))
  }

  if (sizeByIntensity==TRUE) {
    volcanoPlot = volcanoPlot + aes(size=LogInt) +
      scale_size_continuous()
  }

  if (!is.null(pthresholdLine)){
    volcanoPlot = volcanoPlot +
      geom_hline(yintercept=-log10(pthreshold), color=pthresholdLine,
                 alpha=0.5, size=refLineThickness)
  }

  if (!is.null(foldChangeLines)){
    volcanoPlot = volcanoPlot +
      geom_vline(xintercept=foldChangeLines, color=symbolFill["Increased"],
                 alpha=0.5, size=refLineThickness) +
      geom_vline(xintercept=-foldChangeLines, color=symbolFill["Decreased"],
                 alpha=0.5, size=refLineThickness)
  }

  ### Add Labels

  if (is.null(xlab)){ #use colname unless supplied as argument
    volcanoPlot = volcanoPlot + xlab(xlabel)
  } else {
    volcanoPlot = volcanoPlot + xlab(xlab)
  }
  if (is.null(ylab)){
    volcanoPlot = volcanoPlot + ylab(ylabel)
  } else {
    volcanoPlot = volcanoPlot + ylab(ylab)
  }
  if (!is.null(title)){
    volcanoPlot = volcanoPlot +
      ggtitle(title)
  }

  #Set the font size before placing the legend
  if (tolower(themeStyle) == "bw") {
    volcanoPlot = volcanoPlot + bwTheme(baseFontSize) #+ baseTheme(baseFontSize)
  } else {
    volcanoPlot = volcanoPlot + greyTheme(baseFontSize) #+ baseTheme(baseFontSize)
  }

  volcanoPlot = setLegendPosition(volcanoPlot, legendPosition, themeStyle)

  return(volcanoPlot)
}

