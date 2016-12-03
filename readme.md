#DGE.Tools: RNA-Seq Analysis Package

17Aug2016: A new production version (1.0.8) will be released soon.  In the meantime, you should use the dev version 1.0.7.4 because of many
new features and a few bug fixes (in the conversion functions).

DGE.Tools is a suite of functions to facilitate and standardize RNA-Seq DGE analysis.

**Note:** as of Version 1.0 (29Dec2015) to use DGE.Tools you must upgrade to R v3.2.3 (released 10Dec2015) to higher and  ggplot2 v2.0.0 or higher.

v1.0.6 Updates:

* Function runContrasts: Removed the report of genes meeting specified thresholds.  The associated SigTable and SigTablePNG arguments are thus deprecated.  Code using those arguments will still work but those arguments are ignored.  The signature count functionality is now implemented in a separate function.
* Function runContrasts: Qvalue argument is deprecated and ignored.  Qvalue calculations relegated to a separate function
* Function runQvalue added:  Appends Qvalue and qvalue.lfdr columns to topTable data.frames.
* Function runIHW added: Appends an intensity weighted FDR value to topTable data.frames (ihw.adj.P.value) (see the help for reference)
* plotDispersion: added robust=TRUE argument to be consistent with runVoom
* Function EndemblGene2GeneSym added: Uses annotables package to add gene symbol, biotype and description columns to topTable data.frames.


## DGE workflow:

DGE.Tools has modular functions to conduct DGE analysis from counts to contrasts with facility to select detected genes, normalize data (EdgeR TMM), linear modeling (limma voom and lmfit), and contrast analysis (topTable, topTreat). This process is broken logically into steps so that it is easy to, for example, substitute in a new or customized normalization step and still be able to take advantage of the other pieces of the pipeline. Latest additions to the pipeline tools include support for qualityWeights, duplicateCorrelation and SVA analysis.

## Standard data structures for DGE:

Two specialized data structures are employed and provide a standard that will facilitate data sharing.

Raw count data and associated annotation (genes and sample annotation) and other "assay" data (FPKM, CPM, etc) are encapsulated in a SummarizedExperiment. This SummarizedExperiment object (SE) provides a convenient single data object that holds an entire experiment's worth of data. Conveniently, the top level SummarizedExperiment can be subsetted like a dataframe or matrix. Methods implemented for the SummarizedExperiment take care of subsetting all the contained objects properly. While you can manually load data into an SE object, the use of the DGE.Tools functions to create the SE objects insures that the subelements of the SE object will have consistent names and thus facilitate sharing data between analysts.

Alas, the SE is not suitable for capturing all the downstream analysis results so a second data structure has been introduced. Results of linear modeling are encapsulated in a SubsettableListOfArrays (SLOA). The SLOA contains your Design Matrix, Fit object derived from running voom and lmFit.

## QC Plots:

Several QC plots are availble to monitor the quality of your results. These include:

**EdgeR dispersion plot**   
**voom Mean-Variance plot**   
**plotPvalHist**: Faceted plot of pvalue distributions for each contrast to evaluate quality of your Fit.   
**cdfPlot**: Faceted plot of pvalue distributions for each contrast to evaluate quality of your Fit.   

##Data Exploration Plots:

**profilePlot**: Plot LogIntensity vs. LogRatio from topTable dataframes with highlighting of significantly regulated genes.  
**volcanoPlot**: Plot LogRatio vs. NegLogPvalue from topTable dataframes with highlighting of significantly regulated genes.  
**comparePlot**: Compare LogRatios for two samples showing common or uniquely regulated genes  
**ggplotMDS**: Run Multi-Dimentional Scaling analysis and plot the results  

## Scalable Themes:

The plots described below all use scalable fonts in the themes.  Two themes were created for this purpose and are generally useful for your custom ggplotting as well.  The themes greyTheme and bwTheme are available when DGE.Tools is loaded.  They are equivalent to the default ggplot themes theme\_grey and theme\_bw but use scalable fonts.  Another theme, baseFont can be used to scale fonts without disturbing any other plot customizations.  All the themes are functions that take a single argument to set the base font size for a plot.  Finally, printAndSave can be used to print to knitr with a smaller font and simulataneously save the same plot to a image file with a larger font suitable for PowerPoint presentations.

## Other Documentation

See DGE_Tools Rollout.pptx for information on the workflow tools.   
See Build_RSE_Object.pdf for information about creating a RSE and accessing data in a RSE   
See DGE_Workflow_Tutorial.pdf for a detailed tutorial on running various DGE workflows.   
See DGE.ToolsPlotGallery.pdf for code examples and example plots from the data exploration plots. 

## Test Data and Markdown

Example data with associated .Rmd markdown scripts has been posted here:  
/stf/biohtml/thompj27/DGE.Tools/DGE.Tools_Example

Download the entire folder and RData subfolder.  You'll need to modify the setwd line in the first code block of each Rmd file.   

## Installation

CRAN and Bioconductor package dependencies will be installed automatically.  However, we use a few tools not yet 
deposited in those places.  So there are a few pre-requisites before installing DGE.Tools.

DGE.Tools is dependent on the v3.2.3 version of R released 10Dec2015.  Also, make sure to update
ggplot as well.  Probably a good idea to just update all your out of date packages at this point.

It is best to run the install from a fresh R session before loading any packages.  Loaded packages cannot be updated.

```r
    #if you don't have the devtools package already
    install.packages("devtools") 

    #next line required so that missing Bioconductor packages will install
    source("https://bioconductor.org/biocLite.R")

    #you may also need to set proxy information to install from external github accounts
    httr::set_config(httr::use_proxy(url="http://proxy-server.bms.com:8080"))

    #This one provides Ensembl annotation (currently needed for function Xpress2RSE)
    devtools::install_github("stephenturner/annotables")

    #Installs Ron's zFPKM package
    devtools::install_git("http://biogit.pri.bms.com/thompj27/zFPKM")
	
    #now you should have all the dependencies in place to install DGE.Tools 

    #Install the dev version of DGE.Tools 
    devtools::install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools", branch="dev", repos=BiocInstaller::biocinstallRepos()) 
  
```   

After all that, I've still seen some contexts where not all the dependencies are being installed automatically.  Just use
biocLite or install.packages to install the ones it complains are missing.

Finally, There is still one critical piece of R source code that breaks when encapsulated in a
package. The SubsettableListOfArrays.R file still needs to be manually copied to
your local drive and sourced whenever you are manipulating SubsettableListOfArray
objects (i.e. the output of runEdgeRNorm and runVoom functions). 
[Download the source file from this link.](http://bioinformatics.bms.com/active/biohtml/thompj27/DGE.Tools/SubsettableListOfArrays.R)

Place it somewhere on your computer.  For example, I use "~/R/lib/" for my personal source
library ( ~ usually maps to MyDocuments on Windows).  Then you need to source the file before working with a SLOA object.

```r
source("~/R/lib/SubsettableListOfArrays.R")
```
