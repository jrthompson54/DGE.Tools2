### Function printAndSave ###
#' Function  printAndSave
#'
#' Print a ggplot2 object to the console/knitr report and save as a graphic file.
#' Defaults to larger fonts for the graphic file.  Graphic file type is taken from
#' the given filename and can be one of (png, bmp, tiff, jpg, pdf, svg).
#'
#' I found I didn't like the default Legend sizes. I inserted code to increase the
#' size of the print/knitr legend and decrease the size of the legend in the saved
#' graphic file.  This works well for images of about 5-7 inches but may not scale
#' well to other image sizes.  Set scaleLegend to FALSE to get the default behavior
#' back or use legend.key.size in a theme to tweak to your own specifications.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2, png, bmp, tiff, jpeg, pdf
#'
#' @param plotObject A ggplot2 plotobject
#' @param filename A path/filename for the graphic
#' @param width Graphic width in inches (default = 7)
#' @param height Graphic height in inches (default = 5)
#' @param units Units for height and width ("in"|"cm"|"mm") (Default = "in")
#' @param scale Multiplicative scaling factor (Default = 1)
#' @param res Resolution in ppi (default=300)
#' @param printFontSize Base font size for the graphic on the console/knitr (default=12)
#' @param saveFontSize Base font size for the graphic file (default=24)
#' @param scaleLegend Scale the legend smaller if font > 14  (Default = TRUE)
#' @param printPlot Print to console if TRUE (Default=TRUE)
#' @param savePlot Print to file if TRUE (Default = TRUE)
#'
#' @return The print object
#'
#' @examples
#' printAndSave (Myggplot, "myfile.png") #all defaults
#' printAndSave (Myggplot, "myfile.png", width=5, height=4, res=150, printFontSize=10,
#'                saveFontSize=18)  #set a few options
#'
#' @import ggplot2
#' @importFrom grDevices dev.cur dev.off
#' @importFrom assertthat assert_that
#'
#' @export
printAndSave <- function (plotObject, filename, width=7, height=5,
                     units='in', res=300, scale=1,
                     printFontSize=12, saveFontSize=24,
                     scaleLegend = TRUE, printPlot=TRUE, savePlot=TRUE){

  #Save the starting dev level
  startDev <- grDevices::dev.cur()

  #scale the legend text for saved graphics
  save.Legend.ScaledSize <- 10/printFontSize
  save.Legend.ScaledFont <- 7/printFontSize
  print.Legend.ScaledSize <- 0.9 #theme_grey defaults
  print.Legend.ScaledFont <- 0.7

	LegendPrint = theme(
	  legend.text = element_text(colour="Black", size=rel(print.Legend.ScaledFont)),
	  legend.title = element_text(colour="Black",
	                              size=rel((print.Legend.ScaledFont*1.2))),
	  legend.key.size = unit(print.Legend.ScaledSize, "lines"),
	  legend.title.align = 0.5
	)

	LegendSave = theme(
	  legend.text = element_text(colour="Black", size=rel(save.Legend.ScaledFont)),
	  legend.title = element_text(colour="Black",
	                              size=rel((save.Legend.ScaledFont*1.2))),
	  legend.key.size = unit(save.Legend.ScaledSize, "lines"),
	  legend.title.align = 0.5
	)

  #create the output directory if necessary
  if (!file.exists(dirname(filename))) {
    dir.create(dirname(filename), recursive=TRUE)
  }

  #save to file
  if (savePlot==TRUE) {
    #get the file extension
    filetype = tolower(tools::file_ext(filename))

    plot = plotObject

    if (scaleLegend == TRUE && saveFontSize > 14) {
      plot = plot + LegendSave
    }
    supportedFiletypes <- c("png", "bmp", "tiff", "jpeg", "pdf", "svg", "wmf")
    if (filetype %in% supportedFiletypes) {
    	ggsave(filename=basename(filename), plot=plot,
    	       device=filetype, path=dirname(filename),
    	       width=width, height=height, units=units, dpi=res)
    } else {
      warning("Warning: File extension not recognized. No file saved.")
    }
  }

  if (printPlot == TRUE) {
    if (scaleLegend == TRUE) plotObject = plotObject + LegendPrint
    return(plotObject)  #print to console or knitr report
  } else {
    return(NULL)
  }

  #Reset to starting dev level before exit (traps for a malformed ggplot that opens a device and never closes)
  while (grDevices::dev.cur() > startDev)
      grDevices::dev.off()
}

### Function printWithFootnote ###
#' Function  printWithFootnote
#'
#' Print a ggplot2 object to the console/knitr report adding footnote text
#' under the plot.  Use when you want the footnot to be underneath the plot
#' labels.  Only prints the footnote once on a facetted plot.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2, png, bmp, tiff, jpeg, pdf
#'
#' @param plotObject A ggplot2 plotobject
#' @param footnote A path/filename for the graphic
#' @param fontface fontface for the footnotw (default = "plain")
#' @param fontsize size of the footnote font (default = 10)
#' @param hjust  Specify horizontal justification (Default = -0.1)
#'
#' @return Prints the graphic object to the console
#'
#' @examples
#'     #Write to the console or knitr report
#'     printWithFootnote(Myggplot, footnote = "Footnote Text")
#'
#'     #Capture to a file
#'     png ("myplot.png", width=5, height=4, unit="in")
#'     printWithFootnote(Myggplot, footnote = "Footnote Text")
#'     invisible(dev.off())
#'
#' @import ggplot2 tools grDevices
#' @importFrom grid grid.newpage grid.draw textGrob gpar
#' @importFrom assertthat assert_that
#' @importFrom gridExtra arrangeGrob
#'
#' @export
printWithFootnote <- function(plotObject, footnote, fontface="plain", fontsize=10, hjust=-0.1){
  #fontfact values = "plain", "bold", "italic"
  assertthat::assert_that("ggplot" %in% class(plotObject))
  assertthat::assert_that(class(footnote)[[1]] == "character")
  grid::grid.newpage()
  g <- gridExtra::arrangeGrob(plotObject,
                              bottom = grid::textGrob(footnote, x = 0,
                                                      hjust = hjust, vjust=0.1,
                                                      gp = grid::gpar(fontface = fontface, fontsize = fontsize))
                              )
  grid::grid.draw(g)
}
