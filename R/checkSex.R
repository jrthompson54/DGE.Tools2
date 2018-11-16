### Function checkSex ###
#' Function checkSex
#'
#' Take a DGEobj as input and plot expression of XIST vs the highest expressed Y chr gene.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords MDS, RNA-Seq, DGE, QC
#'
#' @param dgeObj A DGEobj (Required)
#' @param species One of "human", "mouse", "rat"
#' @param sexCol Character string name of the sex annotation column in the design table within the dgeObj (optional)
#' @param labelCol Character string name of the design column to use to label points with ggrepel (optional if showLabels=FALSE)
#' @param showLabels Set TRUE to turn on ggrepel labels (default = FALSE)
#' @param chrX Character string name of the X chromosome in the gene annotation within the geneData object (Default = "X")
#' @param chrY Character string name of the Y chromosome in the gene annotation within the geneData object (Default = "Y")
#' @param baseFontSize Sets the base font size for ggplot themes (Default = 14)
#' @param orig Set to FALSE to use your filtered DGEobj data.  Set to TRUE to use original unfiltered data (Default = TRUE)
#'
#' @return A ggplot object
#'
#' @examples
#'
#'    checkSex(dgeObj, species="human")
#'
#' @import ggplot2 magrittr ggrepel ggiraph
#' @importFrom assertthat assert_that
#' @importFrom edgeR calcNormFactors
#' @importFrom dplyr filter select left_join mutate arrange
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @export
checkSex <- function(dgeObj, species, sexCol, labelCol, showLabels=FALSE,
                     chrX = "X", chrY = "Y", baseFontSize=14, orig=TRUE) {

  #convert this to a function in DGE.Tools2
  #input dgeObj, optionally specify x and y genes
  ##if genes not specified, pick highest expressed of each x and y gene
  assertthat::assert_that(!missing(species),
              !missing(dgeObj),
              tolower(species) %in% c("human", "mouse", "rat")
              )

  if (tolower(species) == "rat") { #no XIST in Rat
    x <- .getTopExpressedGene(dgeObj, chr=chrX, orig=TRUE)
  } else {
    x <- switch(tolower(species),
                human = list(gene="ENSG00000229807", genename="XIST"),
                mouse = list(gene="ENSMUSG00000086503", genename="Xist")
    )
  }

  # x <- .getTopExpressedGene(dgeObj, "X")
  y <- .getTopExpressedGene(dgeObj, chr=chrY, orig=TRUE)

  #get data for the two x and y genes
  if (orig){
    log2CPM <- convertCounts(getItem(dgeObj, "counts_orig"), unit="cpm", log=TRUE, normalize="tmm")
  } else {
    log2CPM <- convertCounts(getItem(dgeObj, "counts"), unit="cpm", log=TRUE, normalize="tmm")
  }

  idx <- rownames(log2CPM) %in% c(x$gene, y$gene)
  plotDat <- t(log2CPM[idx,]) %>%
    as.data.frame %>%
    tibble::rownames_to_column(var="rowname")

  #add sample identifier data
  dtemp <- getItem(dgeObj, "design_orig") %>%
    tibble::rownames_to_column(var="rowname")
  #add labelCo and sexCol data as needed
  if (!missing(labelCol) & !missing(sexCol)){
    dtemp %<>% dplyr::select(rowname, labelCol=labelCol, sexCol=sexCol)
    plotDat <- dplyr::left_join(plotDat, dtemp)
    colnames(plotDat) <- c("SampID", x$genename, y$genename, labelCol, sexCol)
    sexPlot <- ggplot(plotDat, aes_string(x=x$genename, y=y$genename, label=labelCol, color=sexCol))
  } else if (!missing(labelCol) & missing(sexCol)){
    dtemp %<>% dplyr::select(rowname, labelCol=labelCol)
    plotDat <- dplyr::left_join(plotDat, dtemp)
    colnames(plotDat) <- c("SampID", x$genename, y$genename, labelCol)
    sexPlot <- ggplot(plotDat, aes_string(x=x$genename, y=y$genename, label=labelCol))
  } else if (missing(labelCol) & !missing(sexCol)){
    dtemp %<>% dplyr::select(rowname, sexCol=sexCol)
    plotDat <- dplyr::left_join(plotDat, dtemp)
    colnames(plotDat) <- c("SampID", x$genename, y$genename, sexCol)
    sexPlot <- ggplot(plotDat, aes_string(x=x$genename, y=y$genename, color=sexCol))
  } else if (missing(labelCol) & missing(sexCol)){
    colnames(plotDat) <- c("SampID", x$genename, y$genename)
    sexPlot <- ggplot(plotDat, aes_string(x=x$genename, y=y$genename))
  }

  sexPlot <- sexPlot +
    geom_point(size=4, shape=21) +
    xlab(paste("X (", x$genename, ")", sep="")) +
    ylab(paste("Y (", y$genename, ")", sep="")) +
    ggtitle("Sex Determination Plot") +
    theme_grey()
  if (showLabels == TRUE & !missing(labelCol))
    sexPlot <-sexPlot + geom_label_repel (label.size=0.125)

  return(sexPlot)
}


#non-exported functions
#
#
#support function for Sex analysis
.getChrDat <- function(dgeObj, chr, orig=FALSE){
  # chr can be one or more chromosome names
  # Returns NULL if none found
  if (orig == TRUE) {
    geneData <- getItem(dgeObj, "geneData_orig")
    log2CPM <- convertCounts(getItem(dgeObj, "counts_orig"), unit="cpm", log=TRUE, normalize="tmm")
  } else {
    geneData <- getItem(dgeObj, "geneData")
    log2CPM <- convertCounts(getItem(dgeObj, "counts"), unit="cpm", log=TRUE, normalize="tmm")
  }

  #filter to specified chr
  if (!is.null(chr))
    geneData %<>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="geneid") %>%
    dplyr::filter(toupper(Chromosome) %in% toupper(chr)) %>%
    tibble::column_to_rownames(var="geneid")

  log2CPM_mat <- log2CPM[rownames(log2CPM) %in% rownames(geneData),]

  if (is.null(nrow(log2CPM_mat))){
    return(NULL)
  } else {
    return(log2CPM_mat)
  }
}

.getTopExpressedGene <- function (dgeObj, chr, orig=FALSE){
  #get top mean expressed gene from given chromosome

  if (orig == TRUE) {
    geneData <- getItem(dgeObj, "geneData_orig")
  } else {
    geneData <- getItem(dgeObj, "geneData")
  }

  #filter geneData to selected chr
  if (!is.null(chr))
    log2CPM_mat <- .getChrDat(dgeObj, chr=chr, orig=orig)
  if (is.null(log2CPM_mat))
    stop ("No gene data found for specified chromosome.")

  meanLogCPM <- rowMeans(log2CPM_mat)

  #get the top expressed gene data
  #Calc mean and sort on mean descending
  log2CPM  <- log2CPM_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="geneid") %>%
    dplyr::mutate(meanLogCPM = meanLogCPM) %>%
    dplyr::arrange(desc(meanLogCPM)) %>%
    dplyr::mutate(meanLogCPM = NULL) %>%
    tibble::column_to_rownames(var="geneid")
  gene <- rownames(log2CPM)[1]
  genename <- geneData$GeneName[rownames(geneData) == gene]
  # dat <- log2CPM %>%
  #   rownames_to_column(var="geneid") %>%
  #   filter(geneid == gene) %>%
  #   column_to_rownames(var="geneid") %>%
  #   as.numeric

  return(list(gene=gene, genename=genename))
}
