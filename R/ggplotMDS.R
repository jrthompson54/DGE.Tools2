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
#' @author John Thompson, \email{jrt@@thompsonclan.org}
#' @keywords MDS, RNA-Seq, DGE, QC
#'
#' @param DGEdata A DGEList object taken after normalization
#'   OR a DGEobj that contains a DGEList OR a log2cpm matrix. (Required)
#' @param colorBy A grouping vector to color by (e.g. ReplicateGroup) (Required)
#' @param shapeBy A grouping vector to map to shape (Optional)
#' @param sizeBy A numeric vector to define point size (Optional)
#' @param top Number of most variant genes to include (default = Inf)
#' @param labels Text labels for the samples.  These should be short
#'   abbreviations of the sample identifiers.
#'   Default = ReplicateGroup or rownames of DGEdata.  Set to NULL to disable
#'   text labels
#' @param labelSize control the Size for the text labels in the plot if you
#'   don't like the default
#' @param title A title for the plot (Optional)
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
#' @param alpha Set transparency (Default = 0.7)
#' @param shapes A vector of shapes to override the default 8 shapes used in shapeBy (optional)
#' @param colors A color pallet to substitute for the default 8 color pallet used by colorBy (optional)
#' @param dim.plot Define which dimension to plot (default = c(1,2))
#' @param shapeName Legend title for shape (optional)
#' @param colorName Legend title for color (optional)
#' @param sizeName Legend title for size (optional)
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
#' @import ggplot2 magrittr ggrepel ggiraph
#' @importFrom assertthat assert_that
#' @importFrom limma plotMDS
#'
#' @export
ggplotMDS <- function(DGEdata,
                      colorBy,
                      shapeBy,
                      sizeBy,
                      top = Inf,
                      labels,
                      labelSize,
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
                      alpha = 0.7,
                      shapes,
                      colors,
                      dim.plot = c(1,2),
                      colorName,
                      shapeName,
                      sizeName
                        )
{


  #default labels to colnames of DGEdata
  addLabels <- TRUE
  if (missing(labels)){ #set default labels
    labels <- colnames(DGEdata)
    # Get labels from ReplicateGroup if present
    if ("DGEobj" %in% class(DGEdata)){
      design <- getItem(DGEdata, "design")
      if (exists("design"))
        if (with (design, exists("ReplicateGroup")))
          labels <- design$ReplicateGroup
    }
  } else if (is.null(labels))
    addLabels <- FALSE

  #argument checks
  if ("DGEobj" %in% class(DGEdata)) #pull out the DGEList
    DGEdata <- getItem(DGEdata, "DGEList")
  else if (!any(c("DGEList", "matrix") %in% class(DGEdata)))
    stop("DGEdata must be class DGEList or DGEobj or matrix")

  assertthat::assert_that(!missing(colorBy),
                          length(colorBy) == ncol(DGEdata))
  if (!missing(shapeBy))
    assertthat::assert_that(length(shapeBy) == ncol(DGEdata))
  if (!missing(sizeBy))
    assertthat::assert_that(length(sizeBy) == ncol(DGEdata))

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
  # if (missing("dim.plot")){
  #   print("dim.plot was missing. dim.plot defaults applied")
  #   dim.plot <- c(1,2)
  # }
  # if (!exists("ndim")){
  #   ndim <- max(dim.plot)
  # }
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

  # grDevices::pdf(NULL) #suppress the plot and just capture the output
  # mds <- plotMDS(DGEdata, top = top, labels = labels, pch = pch,
  #                cex = cex, dim.plot = dim.plot, ndim = ndim,
  #                gene.selection = gene.selection,
  #                xlab = Xlab, ylab = Ylab)
  mds <- limma::plotMDS(DGEdata, top = top, pch = pch,
                 cex = cex, dim.plot = dim.plot, ndim = max(dim.plot),
                 gene.selection = gene.selection,
                 xlab = Xlab, ylab = Ylab,
                 plot=FALSE)
  # invisible(grDevices::dev.off())

  #pull the plotting data together
  if (addLabels == TRUE)
      xydat = data.frame(x=mds$x, y=mds$y, ColorCode=colorBy, Labels=labels)
  else xydat = data.frame(x=mds$x, y=mds$y, ColorCode=colorBy)

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

  #start the plot
  if (byShape == FALSE & bySize == FALSE)
    mdsplot = ggplot(xydat, aes(x=x, y=y, color=ColorCode)) +
        geom_point(shape=symShape, size=symSize, alpha=alpha)

  else if (byShape == TRUE & bySize == FALSE)
    mdsplot = ggplot(xydat, aes(x=x, y=y, color=ColorCode,
                                  shape=Shape)) +
        geom_point(size=symSize, alpha=alpha) +
        scale_shape_manual(values=shapes)

  else if (byShape == FALSE & bySize == TRUE)
    mdsplot = ggplot(xydat, aes(x=x, y=y, color=ColorCode, size=Size)) +
    geom_point(shape=symShape, alpha=alpha)

  else if (byShape == TRUE & bySize == TRUE)
    mdsplot = ggplot(xydat, aes(x=x, y=y, color=ColorCode,
                                shape=Shape,
                                size=Size)) +
    geom_point(alpha=alpha) +
    scale_shape_manual(values=shapes)

  #add point labels
  if (!is.null(labels)){
    if (missing(labelSize))
      mdsplot <- mdsplot +
        geom_text_repel(aes(label = Labels))
    else mdsplot <- mdsplot +
        geom_text_repel(aes(label = Labels), size=labelSize)
  }

  # For discrete color values
  if (length(unique(colorBy)) <= length(colors)){
    mdsplot <- mdsplot +
      scale_fill_manual(values=colors) +
      scale_colour_manual(values=colors)
  }

  #add some other common elements
  mdsplot <- mdsplot +
    coord_fixed() +
    xlab (xylab[[1]]) +
    ylab (xylab[[2]]) +
    ggtitle (title)


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

  #edit legend titles
  if (!missing(colorName)) mdsplot <- mdsplot + labs(color=colorName)
  if (!missing(shapeName) && byShape == TRUE) mdsplot <- mdsplot + labs(shape=shapeName)
  if (!missing(sizeName) && bySize == TRUE) mdsplot <- mdsplot + labs(size=shapeName)

  if (tolower(themeStyle) %in% c("grey", "gray")){
    mdsplot <- mdsplot + theme_grey(baseFontSize)
  } else mdsplot <- mdsplot + theme_bw(baseFontSize)

  return(list(plot=mdsplot, mdsobj=mds))
}

### Function MDS_var_explained ###
#' Function MDS_var_explained
#'
#' Takes a class MDS object from limma::plotMDS and generates two plots: 1)
#' fraction of variance for each dimension, 2) cumulative variance. By default,
#' it plots the first 10 dimensions or the first N dimensions totaling 90%. See
#' topN and cumVarLimit arguments to change the defaults.
#'
#' @author John Thompson, \email{jrt@@thompsonclan.org}
#' @keywords MDS, RNA-Seq, DGE, QC
#'
#' @param mds A class MDS object from limma::plotMDS or a data matrix to analyze
#'   (typically log2) (required)
#' @param topN The number of dimensions to plot (Default = 10)
#' @param cumVarLimit The maximum cumulative variance to plot. Range 0-1. (Default = 0.9)
#' @param barColor Default = "dodgerblue4"
#' @param barFill Default = "dodgerblue3"
#' @param barWidth range 0-1. (Default = 0.65)
#' @param barSize Thickness of the fill border (Default = 0.1)
#' @param baseFontSize Base fontsize for the plot (Default = 14)
#'
#' @return A list with two ggplots and the variance explained data.frame.
#'
#' @examples
#'
#'      #Plot the first two dimensions
#'      MyMDS = ggplotMDS (MyDGEList)
#'      MyMDS[[1]]  #the MDS plot
#'
#'      #Then apply MDS_var_explained (the MDS object is MyMDS[[2]])
#'      varResults <- MDS_var_explained(MyMDS[[2]])
#'      varResults[[1]] #the Variance per dimension plot
#'      varResults[[2]] #the cumulative variance plot
#'      var_explained <- varResults[[3]]  #data used for plotting (unfiltered)
#'
#' @import ggplot2 magrittr
#' @importFrom assertthat assert_that
#' @importFrom limma plotMDS
#' @importFrom matrixStats colVars
#' @importFrom tibble tibble
#'
#' @export
MDS_var_explained <- function(mds,
                              topN = 10,
                              cumVarLimit = 0.9,
                              barColor="dodgerblue4",
                              barFill = "dodgerblue3",
                              barWidth = 0.65,
                              barSize = 0.1,
                              baseFontSize = 14){

  assertthat::assert_that(!missing(mds))

  if (!is(mds, "MDS")) mds <- limma::plotMDS(mds, plot=FALSE)

  mds.distances <- mds %$% distance.matrix %>% as.dist

  mdsvals <- mds.distances %>%
             {suppressWarnings(cmdscale(., k = ncol(mds$distance.matrix)-1))} %>%
             magrittr::set_colnames(stringr::str_c("Dim", seq_len(ncol(.)))) %>%
             as.data.frame

  var_explained <- tibble::tibble(var = mdsvals %>%
                                        as.matrix %>%
                                        matrixStats::colVars() %>%
                                        magrittr::divide_by(sum(.)),
                                  cumvar = cumsum(var),
                                  dim = seq_along(var))

  #often lot's of dimensions.  limit what to plot
  idx <- var_explained$cumvar < cumVarLimit
  if (sum(idx) < topN) topN <- sum(idx)
  plotdat <- var_explained[idx,][1:topN,]

  #plotting
  setBreaks <- function(limits){
    #return integer breaks
    low <- floor(limits[1])
    high <- ceiling(limits[2])
    seq(from=low, to=high, by=1)
  }

  resultList <- list()
  #fractin variance for each dimension
  resultList$varexp <- ggplot(plotdat) +
                       aes(x = dim, y = var) +
                       geom_col(color=barColor,
                                fill=barFill,
                                size=barSize,
                                width=barWidth) +
                       labs(title = "Variance Explained by MDS Dimensions",
                            x = "MDS dimension",
                            y = "Variance Explained") +
                       scale_x_continuous(breaks=setBreaks) +
                       theme_grey(baseFontSize)

  #cumulative variance plot (change the y dimension and relabel)
  resultList$cumvar <- resultList$varexp + aes(y = cumvar) +
                       labs(title = "Cumulative Variance Explained by MDS Dimensions",
                            y = "Cumulative Variance Explained") +
                       ylim(0,1)

  #return the full data table too
  resultList$var_explained <- var_explained

  return(resultList)
}

