#I use this to build and locally install the package while testing.
#Not for use by end Users!

x = getwd()
setwd ("~/R/lib/pkgsrc/DGE.Tools2/")
library(devtools)
document()
load_all()
# build_vignettes()
build()
pkg = "~/R/lib/pkgsrc/DGE.Tools2_0.9.51.tar.gz"
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


library(NCmisc)
flist <- list.functions.in.file("R/runVoom.R", alphabetic=T)
x <- cbind(names(flist), flist) %>% as.data.frame()
x$V1=NULL



> library(DGE.Tools2)
Warning messages:
  1: replacing previous import ‘IRanges::setdiff’ by ‘dplyr::setdiff’ when loading ‘DGEobj’
2: replacing previous import ‘IHW::alpha’ by ‘ggplot2::alpha’ when loading ‘DGE.Tools2’
3: replacing previous import ‘dplyr::combine’ by ‘gridExtra::combine’ when loading ‘DGE.Tools2’
4: replacing previous import ‘assertthat::has_name’ by ‘tibble::has_name’ when loading ‘DGE.Tools2’
5: replacing previous import ‘magrittr::extract’ by ‘tidyr::extract’ when loading ‘DGE.Tools2’




