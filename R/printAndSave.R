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
#' @param res Resolution in ppi (default=300)
#' @param printFontSize Base font size for the graphic on the console/knitr (default=12)
#' @param saveFontSize Base font size for the graphic file (default=24)
#' @param scaleLegend Scale the legend smaller if font > 14  (Default = TRUE)
#' @param printPlot Print to console if TRUE (Default=TRUE)
#' @param savePlot Print to file if TRUE (Default = TRUE)
#'        (Default = TRUE)
#'
#' @return The print object
#'
#' @examples
#' printAndSave (Myggplot, myfile.png) #all defaults
#' printAndSave (Myggplot, myfile.png, width=5, height=4, res=150, printFontSize=10,
#'                saveFontSize=18, style="bw")  #set a few options
#'
#' @import ggplot2 tools grDevices
#'
#' @export
printAndSave <- function (plotObject, filename, width=7, height=5,
                     units='in', res=300, printFontSize=12, saveFontSize=24,
                     scaleLegend = TRUE, printPlot=TRUE, savePlot=TRUE){

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

  if (!file.exists(dirname(filename))) {
    dir.create(dirname(filename), recursive=TRUE)
  }

  filetype = tolower(tools::file_ext(filename))

  if (scaleLegend == TRUE) {
    plotObject = plotObject + LegendPrint
  }

  #save to file
  if (savePlot==TRUE) {
    plot = plotObject

    if (scaleLegend == TRUE && saveFontSize > 14) {
      plot = plot + LegendSave
    }

    if (filetype == "png") {
    	png(file=filename, width=width, height=height, units = units, res = res)
    	print(plot)
    	invisible ( dev.off() )
    } else if (filetype == "bmp") {
      bmp(file=filename, width=width, height=height, units = units, res = res)
      print(plot)
      invisible ( dev.off() )
    } else if (filetype %in% c("tif", "tiff")) {
      tiff(file=filename, width=width, height=height, units = units, res = res, compression="lzw+p")
      print(plot)
      invisible ( dev.off() )
    } else if (filetype %in% c("jpg", "jpeg")) {
      jpeg(file=filename, width=width, height=height, units = units, res = res, quality = 90)
      print(plot)
      invisible ( dev.off() )
    } else if (filetype == "pdf") {
      pdf(file=filename, width=width, height=height, units=units, paper="letter")
      print(plot)
      invisible ( dev.off() )
    } else if (filetype == "svg") {
      svg(file=filename, width=width, height=height, units = units, res = res,
          pointsize = printFontSize)
      print(plot)
      invisible ( dev.off() )
    } else {
      warning("Warning: File extension not recognized. No file saved.")
    }
  }
  
  #make robust against failed plots which leave devices open
  # invisible(
  #     while (dev.cur() >1)
  #       dev.off()
  # )

  if (printPlot == TRUE) {
    return(plotObject)  #print to console or knitr report
  } else {
    return(NULL)
  }
}

