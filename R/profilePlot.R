### Function profilePlot ###
#' Deluxe Profile Plots
#'
#' A profile plot shows Log Intensity on the X axis and LogRatios on the Y axis.
#' This function is intended to show the profile plot from a dataframe
#' created by topTable or topTreat.  Properly normalized data will generally be
#' centered around LogRatio = 0.
#'
#' By default, the plot places "logFC" on the Y axis and "AveExpr" on the X axis.
#' By default, a reference horizontal line is shown at LogRatio = 0 on the Y axis.
#' Optionally, additional reference lines will be drawn at +/- a user supplied
#' logratio threshold. A loess line fit is drawn through the actual data.
#' The points are color coded using the significance measure (pvalue or FDR threshold) supplied by the user.
#' By default, the P.Value field is used with a threshold of 0.01 to color code the points.
#'
#' The profilePlot has also been implemented
#' with a scalable theme that makes it easy to scale font sizes larger
#' or smaller for PPT or Knitr output.  Just add baseFont(n),
#' where n equals the base font size (12 works well for knitr, 18 or 24 works well
#' for PPT).  e.g. MyPlot + theme_grey(18).
#'
#' \strong{Data Structure for the input dataframe:}
#'
#' The defaults are set for dataframes produced by topTable and topTreat.  The columns named "logFC",
#' "AveExpr" and "P.Value" are used by default to accommodate the column
#' names used in topTable/topTreat dataframes.  Any other dataframe
#' can be used with foldchange, intensity and significance measures.  In this case you must use
#' the appropriate arguments to define the column names to use. By default, the
#' column names will be used for the axis labels.  You can overide the xylabels
#' the xlab and ylab arguments.
#'
#' A significance measure (which defaults to P.Value <= 0.01) is used to color code genes that
#' are significantly increased or decreased.  Use the appropriate arguments to use an FDR measure instead.
#'
#' Sensible defaults are chosen for symbols (Size, Shape, Color and Fill).  Optionally,
#' there are arguments to allow you to fiddle with these settings.  A length of 3 is
#' required for these arguments which applies the attributes in this order:
#' Increased, NoChange, Decreased.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords compare plot ggplot2 logratio scatterplot
#'
#' @param df A dataframe with LogRatio and LogIntensity columns and optionally a pvalue or FDR column.
#' @param logRatioCol name of the LogRatio column (Default = "logFC")
#' @param logIntCol name of the LogIntensity column (Default = "AveExpr")
#' @param pvalCol name of the pvalue or FDR column (Default = "P.Value")
#' @param xlab X axis label (default the LogIntensity column name)
#' @param ylab Y axis label (default the LogRatio column name)
#' @param title Plot title (optional)
#' @param pthreshold Used to color points (default = 0.01)
#' @param geneSym Name of the gene symbol column in df.  The gene symbol is
#'    not in topTable output by default so the user has to bind this column
#'    to the dataframe in advance.  Then this column will be used to label
#'    significantly changed points
#' @param rugColor Specify color for rug density plot along xy axes. Set to NULL
#'    to disable rug layer. (Default=NULL)
#' @param rugAlpha Sets the transparency for the rug layer.  An alpha <1.0 takes
#'   a long time to draw (~10 min on 32bit Win7).  So I set the default at 1.0,
#'   but for a final plot try 0.1 or 0.2 for alpha for more informative rug.
#' @param symbolSize Size of symbols for Up, no change and Down. default = c(4, 1, 4);
#'        Note: you cannot choose the exact same size for all three.  But you can use decimal
#'        values.
#' @param symbolShape Shape of the symbols for Up, no change and Down; Default = c(21, 20, 21)
#'        (20 = filled circle, 21 = fillable open circle); Note: you cannot choose the same
#'        symbol shape for all three.
#'        See \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Up, NoChange, Down); default = c("black", "grey25", "grey0")
#'        See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#'        Note: You cannot duplicate colors.
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25 are fillable.
#'        This will have no effect on other symbols.
#'        Default = c("red3", "grey25", "deepskyblue4")
#'        Note: You cannot duplicate colors.
#' @param alpha Controls the transparency of the plotted points (0-1; default = 0.5)
#' @param sizeBySignificance Set to TRUE to size points by the negative Log10 of the
#'        Significance measure (Default = FALSE)
#' @param referenceLine Color for an intercept=0 horizontal reference line
#'        (Default = "blue"; NULL disables)
#' @param foldChangeLines Position of reference horizontal lines for fold change
#'        (Default = log2(1.5); NULL disables)
#' @param lineFitType Enable a linefit through the data (Default = "loess";
#'        "lm" produces a linear fit. NULL disables)
#' @param lineFitColor Color for the fit line (Default = "goldenrod1")
#' @param legendPosition One of "top", "bottom", "left", "right", "ne", "se", "nw", "sw", NULL.
#'        top/bottom/left/right place the legend outside the figure.  ne/se/nw/sw place the figure
#'        inside the figure. NULL disables the legend. Default = "right"
#' @param baseFontSize The smallest size font in the figure in points. Default = 12
#' @param themeStyle "bw" or "grey" which correspond to theme_bw or theme_grey respectively.
#'        Default = bw"
#'
#' @return ggplot object
#'
#' @examples
#'    Simple plot with custom title(df is a topTable dataframe)
#'
#'    MyPlot = profilePlot(df, title = "Plot Title")
#'
#'    Some options with a custom datafile
#'
#'    MyPlot = profilePlot (df, pthreshold = 0.1,
#'      logRatioCol = "Log2ratio",
#'      logIntCol = "AverageIntensity",
#'      pvalCol = "BHFDR",
#'      xlab = "Log2Intensity", ylab = "Log2Ratio",
#'      title = "Profile Plot Title",
#'      referenceLine = "blue",
#'      legendPosition="ne")
#'
#' @import ggplot2 magrittr dplyr
#'
#' @export
profilePlot <- function(df,
                        logRatioCol = "logFC",
                        logIntCol = "AveExpr",
                        pvalCol = "P.Value",
                        pthreshold=0.01,
                        geneSym,
                        rugColor = NULL,
                        rugAlpha = 1.0,
                        xlab=NULL, ylab=NULL, title=NULL,
                        symbolSize = c(4, 1, 4),
                        symbolShape = c(21, 20, 21),
                        symbolColor = c("black", "grey25", "grey0"),
                        symbolFill = c("red3", "grey25", "deepskyblue4"),
                        alpha = 0.5,
                        sizeBySignificance = FALSE,
                        referenceLine = "grey25",
                        foldChangeLines = log2(1.5),
                        refLineThickness = 1,
                        lineFitType = "loess",
                        lineFitColor = "goldenrod1",
                        legendPosition = "right",
                        baseFontSize = 12,
                        themeStyle = "grey"
                        )
{

  #argument checks
  #
  # Make sure specified columns exist
  if (!logRatioCol %in% colnames(df)) {
    stop("LogRatio column not found.")
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

  groupNames <- c("Increased", "No Change", "Decreased")
  names(symbolShape) = groupNames
  names(symbolSize) = groupNames
  names(symbolColor) = groupNames
  names(symbolFill) = groupNames

  ssc = data.frame(group = groupNames,
                   symbolShape = symbolShape,
                   symbolSize = symbolSize,
                   symbolColor = symbolColor,
                   symbolFill = symbolFill,
                   order = c(1,3,2),
                   stringsAsFactors = FALSE) %>% arrange(order)
                    #setting the order defines the legend order

  #columns to plot
  #capture the labels from the colname
  xlabel = logIntCol
  ylabel = logRatioCol
  #now make the columnames suitable for use with aes_string
  x = make.names(colnames(df)[colnames(df)==logIntCol])
  y = make.names(colnames(df)[colnames(df)==logRatioCol])
  colnames(df)[colnames(df)==logIntCol] = make.names(colnames(df)[colnames(df)==logIntCol])
  colnames(df)[colnames(df)==logRatioCol] = make.names(colnames(df)[colnames(df)==logRatioCol])

  #need a NLP column for sizing
  if (sizeBySignificance==TRUE) {
    df$negLog10P = -log10(df[[pvalCol]])
  }

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

  profilePlot <- ggplot (df, aes_string(x=x, y=y)) +
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
      geom_point(alpha=alpha)

###end DELUXE PLOT

  ### Optional Decorations
  if (!is.null(rugColor)){
    profilePlot = profilePlot + geom_rug(data=df, inherit.aes=FALSE,
                                         color = rugColor,
                                         alpha=rugAlpha,
                                         show.legend=FALSE,
                                         aes_string(x=x, y=y))
  }

  if (sizeBySignificance==TRUE) {
    profilePlot = profilePlot + aes(size=negLog10P) +
      scale_size_continuous()
  }

  if (!is.null(referenceLine)){
    profilePlot = profilePlot +
      geom_hline(yintercept=0, color=referenceLine,
                 size=refLineThickness, alpha=0.5)
  }

  if (!is.null(foldChangeLines)){
    profilePlot = profilePlot +
      geom_hline(yintercept=foldChangeLines, color=symbolFill["Increased"],
                 size=refLineThickness, alpha=0.5) +
      geom_hline(yintercept=-foldChangeLines, color=symbolFill["Decreased"],
                 size=refLineThickness, alpha=0.5)
  }

  if (!is.null(lineFitType)){
    profilePlot = profilePlot +
      geom_smooth(aes(group=NULL, shape=NULL, size=NULL, color=NULL, fill=NULL),
                  method = tolower(lineFitType),
                  size=refLineThickness, color=lineFitColor, alpha=alpha,
                  se=FALSE, show.legend=FALSE)
  }

  #Add genesym labels to increased, decreased genes.
  if (!missing(geneSym)){
    #filter df to changed genes
    dfc <- filter(df, group != "No Change")
    profilePlot <- profilePlot +
      geom_text_repel(data=dfc, aes_string(x=x, y=y, label=geneSym),
                      show.legend=FALSE)
  }

  ### Add axis Labels

  if (is.null(xlab)){ #use colname unless supplied as argument
    profilePlot = profilePlot + xlab(xlabel)
  } else {
    profilePlot = profilePlot + xlab(xlab)
  }
  if (is.null(ylab)){
    profilePlot = profilePlot + ylab(ylabel)
  } else {
    profilePlot = profilePlot + ylab(ylab)
  }
  if (!is.null(title)){
    profilePlot = profilePlot +
      ggtitle(title)
  }

  #Set the font size before placing the legend
  if (tolower(themeStyle) == "bw") {
    profilePlot = profilePlot + theme_bw() + baseTheme(baseFontSize)
  } else {
    profilePlot = profilePlot + theme_grey() + baseTheme(baseFontSize)
  }

  profilePlot = setLegendPosition(profilePlot, legendPosition, themeStyle)

  return(profilePlot)
}
