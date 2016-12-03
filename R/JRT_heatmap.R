### Function JRT_heatmap ###
#' Function JRT_heatmap
#'
#' JRT_heatmap is a wrapper around pheatmap.  Pheatmap has about 
#' 40 settable command line arguments.  You might want to set
#' a dozen or so of those options different from the defaults
#' which makes for rather long lines.  But typically, you'll
#' set the same options for several heatmaps.  Enter JRT_heatmap.
#' I redefined the options as a list with sensible defaults
#' established.  For options common to many heatmaps, you only
#' need to set these options once at the top of your script.
#' Then you may need to set no more than the title for each
#' individual heatmap.  
#'
#' Another problem JRT_heatmap solve is uneven scales.  If you
#' scale is assymetric e.g. +6 to -4, then you center color (usually
#' white or black) will not be centered on zero with pheatmap.
#' JRT_heatmap adjusts the scale to insure that the desired center
#' color is centered on zero.
#'
#' Peak at the source code for the full list of settable options.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords heatmap, pheatmap
#'
#' @param mat A numeric matrix of values for the heatmap 
#' @param options A list of pheatmap options 
#'
#' @return a pheatmap plot
#'
#' @examples
#'    #Draw a heatmap with all default values (2D cluster with trees)
#' 		JRT_heatmap (myLogratios)
#'
#'    #Add a title and change a few options
#' 	 	#Here's a title that dynamically adds the number of genes displayed
#'    heatopts$main = paste("My Signature Heatmap\n1D Clusters : ",
#'				nrow(LogRatios), " Genes", sep="")
#'		heatopts$cluster_cols = FALSE
#'		heatopts$clustering_distance_rows = "correlation"
#'		JRT_heatmap(LogRatios, heatopts)
#'
#'		#Save the heatmap to a file instead of printing to the console
#'		heatopts$silent = TRUE
#'		heatopts$filename = "MyHeatmap.png"
#'		JRT_heatmap(LogRatios, heatopts)
#'
#' @import RColorBrewer pheatmap
#'
#' @export
JRT_heatmap <- function(mat, myoptions=NULL){
#JRT 13June2015
# wrapper for pheatmap to simplify calling to:
#    JRT_heatmap(mymatrix, myoptions)
#only need to put non-default options in myoptions list
#
#20June2015
#modify to force zero centered colormap using breaks parameter when breaks is set to zero
#use this when the scale is asymetric to fix the problem that 0 doesnt' map to white or black

# example:  Change two options
# myoptions = list(
#   scale = "column",
#   color = colorRampPalette(c("purple", "blue", "cyan", "green", "yellow", "orange", "red"))(100)
# )
# mymatrix = as.matrix(mtcars)
# JRT_heatmap (mymatrix, myoptions)

#  validoptions = list("color", "kmeans_k", "breaks", "border_color", "cellwidth",
#                      "cellheight", "scale", "cluster_rows", "cluster_cols",
#                      "clustering_distance_rows", "clustering_distance_cols",
#                      "clustering_method", "cutree_rows", "cutree_cols",
#                      "treeheight_row", "treeheight_col", "legend",
#                      "legend_breaks", "legend_labels", "annotation_row",
#                      "annotation_col", "annotation", "annotation_colors",
#                      "annotation_legend", "drop_levels", "show_rownames",
#                      "show_colnames", "main", "fontsize", "fontsize_row",
#                      "fontsize_col", "display_numbers", "number_format",
#                      "number_color", "fontsize_number", "gaps_row",
#                      "gaps_col", "labels_row", "labels_col", "width",
#                      "height", "filename", "silent")
#  library(RColorBrewer)
#  library(pheatmap)
  #Default Options
  options = list(
    color = colorRampPalette(c("blue", "black", "yellow"))(50),
    #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    #colorblind friendly schemes:
    #   Cyan, Black, Red
    #   Cyan, Black, Yellow *JT favorite
    #   Blue, white, yellow
    #   Blue, white, Red  *may project better with dim LCD projectors
    # shades of one color:
    #     colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    kmeans_k = NA,
    breaks = 0,
    border_color = NA, #"grey60",
    cellwidth = NA,
    cellheight = NA,
    scale = "none",  #"row", "column" or "none"
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
        #dist options: "correlation", "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
    clustering_method = "complete",
        #one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
        #"mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
    cutree_rows = NA,
    cutree_cols = NA,
    treeheight_row = 150,
    treeheight_col = 50,
    legend = TRUE,
    legend_breaks = NA,
    legend_labels = NA,
    annotation_row = NA,
    annotation_col = NA,
    annotation = NA,
    annotation_colors = NA,
    annotation_legend = TRUE,
    drop_levels = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = NA,
    fontsize = 10, #points
    fontsize_row = 10, #points
    fontsize_col = 10, #points
    display_numbers = FALSE,
    number_format = "%.2f",  #sprintf style format descriptor
    number_color = "grey30",
    fontsize_number = 8, #points
    gaps_row = NULL,
    gaps_col = NULL,
    labels_row = NULL,  #default is rownames
    labels_col = NULL,  #default is colnames
    width = 10,	#inches
    height = 7,	#inches
    filename = NA,  	#file name for image when silent = T; .png is recommended
    silent = FALSE  #set TRUE to send output to filename
  )

  validoptions = names(options)

  #Swap in the myoptions parameters
  if (!is.null(myoptions)){
    #get the list of custom options
    optnames = names(myoptions)
    #change the default for each option in myoptions
    for (i in 1:length(optnames)){
      if (optnames[i] %in% validoptions){
        options[optnames[i]] = myoptions[optnames[i]]
      } else {
        warning ("option ", optnames[i], " is not a valid option and was ignored!")
      }
    }
  }

  #force zero centered breaks if breaks = 0
  #
  if (!is.na(options$breaks) & length(options$breaks) == 1 & min(mat)<0)  {
    # map even breaks to + and - half
    halflen = round (length(options$color)/2)
#     print(paste("halflen=", halflen))
#     print(paste("min(mat) =", min(mat)))
    upincrement = max(mat)/halflen
    dnincrement = (min(mat)/halflen) %>% abs
#     print(paste("dnincrement =", dnincrement))
    upseq = seq(0+upincrement, max(mat), upincrement)
    dnseq = seq(min(mat), 0-dnincrement, dnincrement)
    options$breaks = c(dnseq, 0, upseq)
#    print (paste("length of Breaks =", length(options$breaks)))
    # length(options$breaks)
  }

  pheatmap (mat,
            color = options$color,
            kmeans_k = options$kmeans_k,
            breaks = options$breaks,
            border_color = options$border_color,
            cellwidth = options$cellwidth,
            cellheight = options$cellheight,
            scale = options$scale,
            cluster_rows = options$cluster_rows,
            cluster_cols = options$cluster_cols,
            clustering_distance_rows = options$clustering_distance_rows,
            clustering_distance_cols = options$clustering_distance_cols,
            clustering_method = options$clustering_method,
            cutree_rows = options$cutree_rows,
            cutree_cols = options$cutree_cols,
            treeheight_row = options$treeheight_row,
            treeheight_col = options$treeheight_col,
            legend = options$legend,
            legend_breaks = options$legend_breaks,
            legend_labels = options$legend_labels,
            annotation_row = options$annotation_row,
            annotation_col = options$annotation_col,
            annotation = options$annotation,
            annotation_colors = options$annotation_colors,
            annotation_legend = options$annotation_legend,
            drop_levels = options$drop_levels,
            show_rownames = options$show_rownames,
            show_colnames = options$show_colnames,
            main = options$main,
            fontsize = options$fontsize,
            fontsize_row = options$fontsize_row,
            fontsize_col = options$fontsize_col,
            display_numbers = options$display_numbers,
            number_format = options$number_format,
            number_color = options$number_color,
            fontsize_number = options$fontsize_number,
            gaps_row = options$gaps_row,
            gaps_col = options$gaps_col,
            labels_row = options$labels_row,
            labels_col = options$labels_col,
            width = options$width,
            height = options$height,
            filename = options$filename,
            silent = options$silent
  )
}
