#DGE.Tools2: RNA-Seq Analysis Package

DGE.Tools2 is a suite of functions to facilitate and standardize RNA-Seq DGE analysis.  DGE.Tools relies on the DGEobj data structure to store DGE data and analysis results.  

**Note:** v0.9.* is the beta version.  All functionality seems to be working but testing is still in progress.  Please let me know
if something isn't working as advertised.  


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

##Data Exploration Plots:

**profilePlot**: Plot LogIntensity vs. LogRatio from topTable dataframes with highlighting of significantly regulated genes.  
**volcanoPlot**: Plot LogRatio vs. NegLogPvalue from topTable dataframes with highlighting of significantly regulated genes.  
**comparePlot**: Compare LogRatios for two samples showing common or uniquely regulated genes  
**ggplotMDS**: Run Multi-Dimentional Scaling analysis and plot the results  


## Other Documentation

See DGE_Tools Rollout 2017.pptx for information on the workflow tools.   
  
See DGE.ToolsPlotGallery.pdf for code examples and example plots from the data exploration plots. 

## Test Data and Markdown

A test dataset can be found in the installed library folder (../library/DGE.Tools2/extdata).
To load the test data:

```
rawDataPath <- paste(.libPaths()[[1]], "/DGE.Tools2/extdata", sep="")
dgeObj <- OmicsoftToDgeObj(path=rawDataPath)
```

A sample analysis markdown is being prepared using this dataset and will be
added to the package as a vignette.  The draft (and fairly complete) version
is [found here](https://biogit.pri.bms.com/thompj27/DGE.Tools2/blob/master/vignettes/DGE.Tools2_TestWorkflow.Rmd)

## Installation (New!)

It is best to run the install from a fresh R session before loading any
packages because loaded packages cannot be updated.

Charles Tilford has created a BMS repository.  Packages from the BMS repository
can now be installed simply with install.packages (after a little setup).

One time setup :

```
    #Install Charles' bmsPackageTools to enable the BMS BRAN 
    #repository (modifies your ~/.Rprofile file)
    install.packages("bmsPackageTools", repos="http://kraken.pri.bms.com/biohtml/nositeminder/BRAN/")
    library('bmsPackageTools')
    setBmsRepo(rprofile=TRUE)

    #if you don't have the devtools package already
    install.packages("devtools")
```

CRAN and Bioconductor package dependencies should be installed automatically. 
However, we use a few tools not yet deposited in those places.  So there are a
few pre-requisites before installing DGE.Tools.

```
    #GIT-based packages

    # You may need to set proxy information to install from external github accounts
    httr::set_config(httr::use_proxy(url="http://proxy-server.bms.com:8080"))
    
    #This one provides Ensembl annotation (used by a few gene ID conversion functions)  
    devtools::install_github("stephenturner/annotables")
```

Then install or update DGE.Tools and related packages with the following commands:  

```
    install.packages("JRTutil")
    install.packages("zFPKM")
    install.packages("DGEobj")
    install.packages("DGE.Tools2")
    install.packages("Xpress2R")
```

# Installation directly from BioGit (Backup method)

You can also install directly from BioGit if the BMS repository is down for some
reason.

```
    #next line required so that missing Bioconductor packages will install  
    source("https://bioconductor.org/biocLite.R")

    #you may also need to set proxy information to install from external github accounts  
    httr::set_config(httr::use_proxy(url="http://proxy-server.bms.com:8080"))

    #This one provides Ensembl annotation (used by a few gene ID conversion functions)  
    devtools::install_github("stephenturner/annotables")  

    #Install Xpress2R package if you wish to retrieve Xpress data as a DGEobj  
    devtools::install_git("http://biogit.pri.bms.com/thompj27/Xpress2R",
                           repos=BiocInstaller::biocinstallRepos())   

    #Install the zFPKM package
    devtools::install_git("http://biogit.pri.bms.com/thompj27/zFPKM")
	
    #now you should have all the dependencies in place to install DGE.Tools 

    #Install the DGEobj Package (The S3 class data structure used by DGE.Tools2
    devtools::install_git("http://biogit.pri.bms.com/thompj27/DGEobj",
                            repos=BiocInstaller::biocinstallRepos())   

    #Install the DGE.Tools2
    devtools::install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools2",
                            repos=BiocInstaller::biocinstallRepos()) 
```   


