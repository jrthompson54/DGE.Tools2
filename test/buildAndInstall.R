#I use this to build and locally install the package while testing.
#Not intended for use by end Users!

x = getwd()
setwd ("~/R/lib/pkgsrc/DGE.Tools2/")
library(devtools)
document()
load_all()
# build_vignettes()
build()
pkg = "~/R/lib/pkgsrc/DGE.Tools2_0.9.69.tar.gz"
install.packages(pkg, repos=NULL, type="source")
setwd(x)

# #install from Git
# #After pushing to git...
# install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools2")
# install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools2", repos=BiocInstaller::biocinstallRepos())
#
# #for dev
# install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools", branch="dev", repos=BiocInstaller::biocinstallRepos())
#


# list all function calls in the input R script.  Useful to search your code for
# function calls that may need to be double colon referenced.
library(NCmisc)
flist <- list.functions.in.file("R/runVoom.R", alphabetic=T)
x <- cbind(names(flist), flist) %>% as.data.frame()
x$V1=NULL


