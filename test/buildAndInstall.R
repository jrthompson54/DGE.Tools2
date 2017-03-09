#I use this to build and locally install the package while testing.
#Not for use by end Users!

x = getwd()
setwd ("~/R/lib/pkgsrc/DGE.Tools2/")
library(devtools)
document()
load_all()
#build_vignettes()
build()
pkg = "~/R/lib/pkgsrc/DGE.Tools2_0.9.1.tar.gz"
install.packages(pkg, repos=NULL, type="source")
setwd(x)

#install from Git
#After pushing to git...
install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools2")
install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools2", repos=BiocInstaller::biocinstallRepos())

#for dev
install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools", branch="dev", repos=BiocInstaller::biocinstallRepos())

#DGEobj
x = getwd()
setwd ("~/R/lib/pkgsrc/DGEobj/")
library(devtools);document();load_all;build()
pkg = "~/R/lib/pkgsrc/DGEobj_0.2.0.tar.gz"
install.packages(pkg, repos=NULL, type="source")
setwd(x)


function todo

	runQvalue.R
	runIHW.R

	JRT_heatmap.R (rewrite to use a list of parameters)
	extractCol.R  (rewrite to check for matching rownames/order)
	
	subset not working with logical idx

INC1585521  Avaya softephone ticket

Next:  IHW and qvalue handling  
These functions work on a contrast list.  
provide a dgeObj and a contrast list (named)
add IHW and qvalue as attributes and replace/overwrite the existing contrast objects.

Stopping point:
    testing runQvalue, addItems
But 
> dgeObj <- runEdgeRNorm(dgeObj)
Error in .isAllZero(counts) : 
    count matrix must be integer or double-precision

