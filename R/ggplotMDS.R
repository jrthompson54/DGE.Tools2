### Function ggplotMDS ###
#' Function ggplotMDS
#'
#' This is a wrapper around the plotMDS function that generates the plot with
#' ggplot2 instead of base graphics.
#'
#' colorBy, shapeBy and sizeBy are grouping variables to encode group info by
#' color, shape or size.  These are vectors that must be the same length as 
#' ncol(DGEdata).  colorBy and sizeBy will plot as continuous color or size
#' changes if a numeric vector is used.  Convert the vector to a factor to 
#' treat as groups instead of continuous.
#' 
#' The underlying plotMDS function uses a default of top=500 to use the top 500
#' highest fold change genes for the analysis.  The method runs quick and I've
#' found situations where a larger value produced a more stable result.  Thus
#' I've change the default to top=Inf.  You should play with this parameter and
#' at least try a number that is close to the number of differential genes
#' detected in your data.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords MDS, RNA-Seq, DGE, QC
#'
#' @param DGEdata A DGEList object taken after normalization 
#'   or a DGEobj that contains a DGEList. (Required)
#' @param colorBy A grouping vector to color by (e.g. ReplicateGroup) (Required)
#' @param shapeBy A grouping vector to map to shape (Optional)
#' @param sizeBy A numeric vector to define point size (Optional)
#' @param top Number of most variant genes to include (default = Inf)
#' @param labels Text labels for the samples.  These should be short
#'   abbreviations of the sample identifiers so that the plot is not too
#'   cluttered.  
#' @param labelSize Size for the text labels in the plot (default = 2)
#' @param textColor Color for the text labels in the plot (default = "blue2")
#' @param vlineIntercept X intercept of vertical line (Optional)
#' @param hlineIntercept Y intercept of horizontal line (Optional)
#' @param reflineColor Color for the horizontal and vertiacl reference lines
#'   (Default = "darkgoldenrod1")
#' @param reflineSize Thickness of the reference lines (Default = 0.5)
#' @param baseFontSize Base fontsize for the plot (Default = 12)
#'
#' @param themeStyle One of "grey" or "bw" (Default = "grey")
#' @param symShape Set the default shape of the symbols if not mapped to a column (Default = 19 solid circle)
#' @param symSize Set the default size of the symbols if not mapped to a column 
#'   (Default = 5)
#' @param symFill Set color for the fill on open symbols (Default = "blue2")
#' @param symColor set color for solid symbols or outline for open symbols
#'   (Default = "blue2")
#' @shapes A vector of shapes to override the default 8 shapes used in shapeBy (optional)
#' @colors A color pallet to substitute for the default 8 color pallet used by colorBy (optional)
#' @param ... arguments passed through to plotMDS function (see ?plotMDS in package
#'    edgeR)
#'
#' @return A list with two elements, the ggplot object and the mds object returned
#'    by the plotMDS function.
#'
#' @examples
#'
#'      #Plot the first two dimensions using all genes
#'      MyMDS = ggplotMDS (MyDGEList)
#'      #display the plot and save a png file.
#'      printAndSave (MyMDS[[1]], file = "MyMDS.PNG", width=5, height=5)
#'
#'      #plot the 2nd and 3rd dimensions using the top 1000 genes
#'      MyMDS = ggplotMDS (MyDGEList, dim.plot=c(2,3) ndim =3)
#'      MyMDS[[1]]
#'
#' @import ggplot2 magrittr edgeR assertthat
#'
#' @export
ggplotMDS <- function(DGEdata,
                      colorBy, 
                      shapeBy,
                      sizeBy,
                      top = Inf,
                      labels,
                      labelSize = 2,
                      title,
                      textColor = "blue2",
                      hlineIntercept,
                      vlineIntercept,
                      reflineColor = "darkgoldenrod1",
                      reflineSize = 0.5,
                      baseFontSize = 12,
                      themeStyle = "grey",
                      symShape = 16, #solid circle
                      symSize = 5,
                      symFill = "blue2",
                      symColor = "blue2",
                      shapes,
                      colors,
                      ...
                        )
{
    #logic 
    #if labels and sym arguments set: plot points with labels
    #if only sym arguments: plot points
    #if only labels: plot text only
  
    textMode <- FALSE
    if (!missing(labels))
        textMode <- TRUE
              
  #argument checks
  if (class(DGEdata)[[1]] == "DGEobj") #pull out the DGEList
      DGEdata <- getItem(DGEdata, "DGEList") 
  else if (class(DGEdata)[[1]] != "DGEList")
    stop("DGEdata must be class DGEList or DGEobj")
  
  assert_that(!missing(colorBy),
              length(colorBy) == ncol(DGEdata))
  if (!missing(shapeBy))
      assert_that(length(shapeBy) == ncol(DGEdata))
  if (!missing(sizeBy))
      assert_that(length(sizeBy) == ncol(DGEdata))
  
  #shapes: solid circle, square, triangle, diamond, open circle, square, triangle, diamond
  myShapes = c(16, 15, 17, 18, 21, 22, 24, 23)
  if (missing(shapes))
      shapes <- myShapes
  
  # ColorBlind palette: 
  # http://www.ucl.ac.uk/~zctpep9/Archived%20webpages/Cookbook%20for%20R%20%C2%BB%20Colors%20(ggplot2).htm
  # 
  cbbPalette <- c("#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#E69F00",  "#F0E442", "#000000")
  if (missing(colors))
      colors <- cbbPalette
    
  #defaults for plotMDS arguments
  if (!exists("dim.plot")){
    dim.plot <- c(1,2)
  }
  if (!exists("ndim")){
    ndim <- 2
  }
  if (!exists("gene.selection")){
    gene.selection <- "pairwise"
  }
  if (!exists("Xlab")){
    Xlab <- NULL
  }
  if (!exists("Ylab")){
    Ylab <- NULL
  }
  if (!exists("pch")){
    pch <- NULL
  }
  if (!exists("cex")){
    cex <- 1
  }
  if(!exists("method")){
    method <- "logFC"
  }
  if (!exists("prior.count")){
    prior.count <- 2
  }
  if (missing(labels)){
    labels <- colnames(DGEdata)
  }
  if (missing(title)){
    title <- "MDS Plot"
  }

  pdf(NULL) #suppress the plot and just capture the output
  mds <- plotMDS(DGEdata, top = top, labels = labels, pch = pch,
                 cex = cex, dim.plot = dim.plot, ndim = ndim,
                 gene.selection = gene.selection,
                 xlab = Xlab, ylab = Ylab)
  invisible(dev.off())

  #pull the plotting data together
  xydat = data.frame(x=mds$x, y=mds$y, ColorCode=colorBy) 
  
  byShape <- FALSE
  if (!missing(shapeBy)){
      xydat$Shape <- shapeBy
      byShape <- TRUE
  }
  bySize <- FALSE
  if (!missing(sizeBy)){
      xydat$Size <- sizeBy
      bySize <- TRUE
  }
  
  xylab = list(paste (mds$axislabel, mds$dim.plot[[1]], sep=" "),
               paste (mds$axislabel, mds$dim.plot[[2]], sep=" "))
  if (!is.null(Xlab))
      xylab[[1]] <- Xlab
  if (!is.null(Ylab))
      xylab[[2]] <- Ylab
  
  my_gg <- g + geom_point_interactive(aes(tooltip = model), size = 2) 
  ggiraph(code = print(my_gg), width = .7)
  
  #start the plot
  if (byShape == FALSE & bySize == FALSE)
    mdsplot = ggplot(xydat, aes(x=x, y=y, color=ColorCode)) +
        geom_point(shape=symShape, size=symSize) +
        scale_fill_manual(values=colors) +
        scale_colour_manual(values=colors) 
  
  else if (byShape == TRUE & bySize == FALSE)  
    mdsplot = ggplot(xydat, aes(x=x, y=y, color=ColorCode, 
                                  shape=Shape)) +
        geom_point(size=symSize) +
        scale_shape_manual(values=shapes) +
        scale_fill_manual(values=colors) +
        scale_colour_manual(values=colors) 
  
  else if (byShape == FALSE & bySize == TRUE)  
    mdsplot = ggplot(xydat, aes(x=x, y=y, color=ColorCode, 
                                  size=Size)) +
        geom_point(shape=symShape) +
        scale_fill_manual(values=colors) +
        scale_colour_manual(values=colors) 
  
  else if (byShape == TRUE & bySize == TRUE)  
    mdsplot = ggplot(xydat, aes(x=x, y=y, color=ColorCode, 
                                  shape=Shape,
                                  size=Size)) +
        geom_point() +
        scale_shape_manual(values=shapes) +
        scale_fill_manual(values=colors) +
        scale_colour_manual(values=colors) 
 
  #add some labels
  mdsplot <- mdsplot + 
      xlab (xylab[[1]]) +
      ylab (xylab[[2]]) +
      ggtitle (title)
  
  
  browser()
  
  #place an annotation on the bottom left of the plot
  xrange <- getXrange(mdsplot)
  yrange <- getYrange(mdsplot)
  #put the annotation 10% from xmin
  xpos <- xrange[1] + ((xrange[2] - xrange[1]) * 0.1 )
  alabel <- paste("top ", mds$top, " genes : gene.selection = ",
                 mds$gene.selection, sep="")
  mdsplot <- mdsplot + annotate ("text", x = xpos, y = yrange[1],
                                label = alabel, hjust=0,
                                size=rel(2.5), color="grey30")

  if (!missing(hlineIntercept)){
    mdsplot <- mdsplot + geom_hline (yintercept = hlineIntercept,
                                     color = reflineColor,
                                     size=reflineSize)
  }
  if (!missing(vlineIntercept)){
    mdsplot <- mdsplot + geom_vline (xintercept = vlineIntercept,
                                     color = reflineColor,
                                     size=reflineSize)
  }

  if (tolower(themeStyle) %in% c("grey", "gray")){
    mdsplot <- mdsplot + theme_grey(baseFontSize)
  } else mdsplot <- mdsplot + theme_bw(baseFontSize)

  return(list(plot=mdsplot, mdsobj=mds))

}

