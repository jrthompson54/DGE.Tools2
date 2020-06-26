# DGE.Tools2: RNA-Seq Analysis Workflow Package

DGE.Tools2 is a suite of functions to facilitate and standardize RNA-Seq DGE analysis.  DGE.Tools relies on the DGEobj data structure to store DGE data and analysis results.  


## DGE workflow:

DGE.Tools2 has modular functions to conduct DGE analysis from counts to
contrasts with facility to select detected genes, normalize data (EdgeR TMM),
linear modeling (limma voom and lmFit), and contrast analysis (topTable,
topTreat). This process is broken logically into steps so that it is easy to,
for example, substitute in a new or customized normalization step and still be
able to take advantage of the other pieces of the pipeline. The DGE.Tools
workflow include support for qualityWeights, duplicateCorrelation and SVA
analysis.  The most important reason  to use DGE.Tools2 is that it produces a
standardized data object, the DGEobj, that captures and annotates your
workstream making your data better documented and easier to incorporate into
downstream integrative analyses.

## Standard data structures for DGE:

DGE.Tools2 uses an S3 data class called DGEobj to capture results in a customizable and reusable data object.  
See the [DGEobj package](https://biogit.pri.bms.com/thompj27/DGE.Tools2) for more details.

## QC Plots:

Several QC plots are availble to monitor the quality of your results. These include:

**EdgeR dispersion plot**   
**voom Mean-Variance plot**   
**plotPvalHist**: Faceted plot of pvalue distributions for each contrast to evaluate quality of your Fit.   
**cdfPlot**: Faceted plot of pvalue distributions for each contrast to evaluate quality of your Fit. 
**QCplots**: Plot alignment metrics from Omicsoft or Xpress  

## Data Exploration Plots:

**profilePlot**: Plot LogIntensity vs. LogRatio from topTable dataframes with highlighting of significantly regulated genes.  
**volcanoPlot**: Plot LogRatio vs. NegLogPvalue from topTable dataframes with highlighting of significantly regulated genes.  
**comparePlot**: Compare LogRatios for two samples showing common or uniquely regulated genes  
**ggplotMDS**: Run Multi-Dimentional Scaling analysis and plot the results  


## Other Documentation

* DGE_Tools_Training_Mar2019.pptx   
* vignettes/DGE.Tools2_Workflow.pdf:  Workflow example  
* vignettes/DGE.ToolsPlotGallery.pdf: code examples for data exploration plots  


## Installation 

It is best to run the install from a fresh R session before loading any
packages because loaded packages cannot be updated.

 

```r
    require(devtools)
    devtools::install_git("https://github.com/jrthompson54/DGE.Tools2") 
```   
