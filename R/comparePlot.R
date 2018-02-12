### Function comparePlot ###
#' Deluxe Compare Plots
#'
#' In the simpliest form it draws a nicely formatted scatterplot of the
#' first two columns in a dataframe.  Several formatting options are
#' available to enhance the plots.  If you supply pvalues or FDR values and a threshold,
#' the plot is color coded to show X unique, Y Unique and
#' Common differentialy expressed (DE) genes in different colors.
#'
#' Other options add
#' an identity line (slope=1, intercept=0) and/or a crosshair at (0,0).
#'
#' The comparePlot has also been implemented
#' with a scalable theme that makes it easy to scale font sizes larger
#' or smaller for PPT or Knitr output.  Just add either theme_grey(n) or theme_bw(n),
#' where n equals the base font size (12 works well for knitr, 18 or 24 works well
#' for PPT).  e.g. MyPlot + theme_grey(18).  However, theme_gray and theme_bw have the
#' side effect of resetting the legend position.  Thus, do not use those standard themes
#' if you want to preserve the custom legend locations.  You must rely on the baseFontSize
#' argument in this function to set font sizes to preserve a custom legend location.
#'
#' \strong{Data Structure for the input dataframe:}
#'
#' The x and y values should be in the first two columns. By default, their
#' column names will be used for the axis labels.  You can overide the xylabels
#' the xlab and ylab arguments.
#'
#' Optionally, significance measures in the form of pvalues or FDR values can be supplied
#' for X and Y respectively.  If provided, these columns \strong{must} be named "xp" and "yp" respectively.
#' Together with a threshold (which defaults to 0.01), these values will
#' be used to color code the plot to show X Unique, Y Unique and Common DE
#' genes.  You can use either pvalues or FDR values for the signficance columns.  Use the
#' pthreshold argument to set a proper threshold for FDR values.
#'
#' Sensible defaults are chosen for symbols (Size, Shape, Color and Fill).  Optionally,
#' there are arguments to allow you to fiddle with these settings.  A length of 4 is
#' required for these arguments which applies the attributes in this order:
#' Common, X Unique, Y Unique and Not Significant.
#'
#' Note if you're not using pvalues or FDR values to color the plot, the X Unique color
#' values are used.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords compare plot ggplot2 logratio scatterplot
#'
#' @param df A dataframe with the first two columns representing the x and y variables.
#'          Optionally add xp and yp columns to hold pvalues or FDR values.
#' @param xlab X axis label (default to first column name)
#' @param ylab Y axis label (default to second column name)
#' @param title Plot title (optional)
#' @param pthreshold Used to color points (default = 0.01)
#' @param symbolSize Size of symbols; default = c(4, 4, 4, 2)
#' @param symbolShape Shape of the symbols; Default = c(21, 21, 21, 20)  (20 = filled circle,
#'        21 = fillable open circle) See \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Common, xUnique, yUnique, NoChange);
#'        default = c("black", "grey0", "grey1", "grey25")
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25 are fillable. This will
#'        have no effect on other symbols. Default = c("darkgoldenrod1", "deepskyblue4", "red3", "grey25")
#' @param alpha Controls the transparency of the plotted points (0-1; default = 0.5)
#' @param rugColor Specify color for rug density plot along xy axes. Set to NULL
#'    to disable rug layer. (Default=NULL)
#' @param rugAlpha Sets the transparency for the rug layer.  An alpha <1.0 takes
#'   a long time to draw (~10 min on 32bit Win7).  So I set the default at 1.0,
#'   but for a final plot try 0.1 or 0.2 for alpha for more informative rug.
#' @param crosshair Color for the crosshair; (Default = "grey50", NULL disables)
#'        See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#' @param referenceLine Color for a slope=1, intercept=0 reference line
#'        (Default = "darkgoldenrod1"; NULL disables)
#' @param refLineThickness Set thickness for crosshair and referenceLine (Default=1)
#' @param dens2D Overlay a density2D layer (Default = TRUE) Note: only applies
#'        to simple compare plot without pvalue coloring.
#' @param legendPosition One of "top", "bottom", "left", "right", "ne", "se", "nw", "sw", NULL.
#'        top/bottom/left/right place the legend outside the figure.  ne/se/nw/sw place the figure
#'        inside the figure. NULL disables the legend. Default = "right"
#' @param baseFontSize The smallest size font in the figure in points. Default = 12
#' @param themeStyle "bw" or "grey" which corresponds to theme_bw or theme_grey respectively.
#'        Default = bw"
#' @param footnote optional string placed right justified at bottom of plot.
#' @param footnoteSize applies to footnote. (default = 3)
#' @param footnoteColor applies to footnote. (default = "black")
#' @param footnoteJust Value 0-1. 0 is left justified, 1 is right justified, 0.5 is centered. (default=1)
#'
#' @return ggplot object
#'
#' @examples
#'    All defaults with a custom title
#'
#'    MyPlot = comparePlot(df, title = "Plot Title")
#'
#'    Deluxe Plot with all the bells and whistles. Note df contains columns xp and yp.
#'
#'    MyPlot = comparePlot (df, pthreshold = 0.5,
#'      xlab = "x Axis Label", ylab = "y Axis Label",
#'      title = "Plot Title",
#'      crosshair = "red",
#'      referenceLine = "blue",
#'      legendPosition="nw")
#'
#' @import ggplot2 magrittr dplyr
#'
#' @export
comparePlot <- function(df, pthreshold=0.01,
                        xlab=NULL, ylab=NULL, title=NULL,
                        symbolSize = c(4, 4, 4, 2),
                        symbolShape = c(21, 21, 21, 20),
                        symbolColor = c("black", "grey0", "grey1", "grey25"),
                        symbolFill = c("darkgoldenrod1", "deepskyblue4", "red3", "grey25"),
                        alpha = 0.5,
                        rugColor = NULL,
                        rugAlpha = 1.0,
                        crosshair="grey50",
                        referenceLine="darkgoldenrod1",
                        refLineThickness = 1,
                        dens2D = TRUE,
                        legendPosition = "right",
                        baseFontSize = 12,
                        themeStyle = "bw",
                        footnote,
                        footnoteSize=3,
                        footnoteColor="black",
                        footnoteJust=1
                        ) {

  with.seed <- function(seed, expr) {
    saved.seed <- .Random.seed
    on.exit(.Random.seed <<- saved.seed)
    set.seed(seed)
    expr
  }

  #argument checks
  #
  # 1st two columns should be numeric
  if (ncol(df) <2) {
    stop ("Need at least 2 numeric columns to proceed.")
  } else if (!is.numeric(df[,1]) || !is.numeric(df[,2])) {
    stop("Column 1 and 2 must hold numeric X and Y values.")
  }
  #symbol parameters must all be length=4
  if (!length(symbolSize)==4 || !length(symbolShape)==4 ||
      !length(symbolColor)==4 || !length(symbolFill)==4){
    stop("symbol arguments must be length=4 (symbolSize, SymbolShape, SymbolColor, symbolFill)")
  }

  names(symbolShape) = c("Common", "xUnique", "yUnique", "NoChange" )
  names(symbolSize) = c("Common", "xUnique", "yUnique", "NoChange" )
  names(symbolColor) = c("Common", "xUnique", "yUnique", "NoChange" )
  names(symbolFill) = c("Common", "xUnique", "yUnique", "NoChange" )

  ssc = data.frame(group = c("Common", "X Unique", "Y Unique", "Not Significant"),
                   symbolShape = symbolShape,
                   symbolSize = symbolSize,
                   symbolColor = symbolColor,
                   symbolFill = symbolFill,
                   stringsAsFactors = FALSE)

  #used to set uniform square scale
  scalemax = df[,1:2] %>% as.matrix %>% abs %>% max %>% multiply_by(1.05)

  #columns to plot
  #capture the labels from the colname
  xlabel = colnames(df)[1]
  ylabel = colnames(df)[2]
  #now make the columnames suitable for use with aes_string
  x = make.names(colnames(df)[1])
  y = make.names(colnames(df)[2])
  colnames(df)[1:2] = make.names(colnames(df)[1:2])

  ###SIMPLE PLOT: plot all data (no significance measures supplied)
  if (is.null(df[["xp"]]) | is.null(df[["yp"]])) {

    CompPlot = ggplot(df, aes_string(x=x, y=y)) +
      geom_point(shape=1,
                 size=symbolSize[["xUnique"]],
                 color=symbolFill[["xUnique"]],
                 fill=symbolFill[["xUnique"]],
                 alpha=alpha) +
      coord_equal(xlim=c(-scalemax, scalemax), ylim=c(-scalemax, scalemax))
    if (dens2D == TRUE){
      CompPlot = CompPlot + geom_density2d(color=symbolFill[["yUnique"]])
    }
  } else
  ### DELUXE PLOT: plot groups in different colors/shapes
  if (!is.null(df[["xp"]]) & !is.null(df[["yp"]])){

    #Let's plot the subsets
    xindx = df[["xp"]] <= pthreshold
    yindx = df[["yp"]] <= pthreshold

    #boolean indexes to parse groups
    bothindx = xindx & yindx
    neitherindx = !xindx & !yindx
    xindx = xindx & !bothindx #unique to X
    yindx = yindx & !bothindx #unique to y

    #create group factor column in df
    df$group=NA
    df$group[bothindx] = "Common"
    df$group[xindx] = "X Unique"
    df$group[yindx] = "Y Unique"
    df$group[neitherindx] = "Not Significant"
    df %<>% left_join(ssc)
    df$group %<>% factor(levels=c("Common", "X Unique", "Y Unique", "Not Significant"))

    #set an order field to control order of plotting
    #plot order is 1) Not Significant plotted first, 2) randomly plot X and Y Unique
    #then plot Common last
    df$order = with.seed(1954, sample.int(nrow(df))) + #add random order (seeded for reproducibility)
      nrow(df) * (df$group == "Common") +  #assign common a high value to sort last
      -nrow(df) * (df$group == "Not Significant") #assign NotSig group neg values to sort first

    CompPlot <- ggplot (df, aes_string(x=x, y=y)) +
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
        #make it square with same axis scales
        coord_equal(xlim=c(-scalemax, scalemax), ylim=c(-scalemax, scalemax)) +
        #box around the legend
        theme(legend.background = element_rect(fill="gray95", size=.5, linetype="dotted"))
  } ###end DELUXE PLOT

  ### Optional Decorations

  if (!is.null(rugColor)){
    CompPlot = CompPlot + geom_rug(data=df, inherit.aes=FALSE,
                                         color = rugColor,
                                         alpha=rugAlpha,
                                         show.legend=FALSE,
                                         aes_string(x=x, y=y))
  }

  if (!is.null(crosshair)){
    CompPlot = CompPlot +
    geom_hline(yintercept=0, color=crosshair,
               size=refLineThickness, alpha=0.5) +
    geom_vline(xintercept=0, color=crosshair,
               size=refLineThickness, alpha=0.5)
  }

  if (!is.null(referenceLine)){
    CompPlot = CompPlot +
      geom_abline(slope=1, intercept=0, color=referenceLine,
                  size=refLineThickness, alpha=0.5)
  }

  ### Add Labels

  if (is.null(xlab)){ #use colname unless supplied as argument
    CompPlot = CompPlot + xlab(xlabel)
  } else {
    CompPlot = CompPlot + xlab(xlab)
  }
  if (is.null(ylab)){
    CompPlot = CompPlot + ylab(ylabel)
  } else {
    CompPlot = CompPlot + ylab(ylab)
  }
  if (!is.null(title)){
    CompPlot = CompPlot +
      ggtitle(title)
  }

  #Set the font size before placing the legend
  if (tolower(themeStyle) == "bw") {
    CompPlot = CompPlot + theme_bw() + baseTheme(baseFontSize)
  } else {
    CompPlot = CompPlot + theme_grey() + baseTheme(baseFontSize)
  }

  CompPlot = setLegendPosition(CompPlot, legendPosition, themeStyle)

  #footnote
  if (!missing(footnote))
    CompPlot <- addFootnote (CompPlot, footnoteText=footnote,
                             footnoteSize=footnoteSize,
                             footnoteColor="black",
                             footnoteJust=footnoteJust,
                             yoffset = 0.05)

  return(CompPlot)
}

