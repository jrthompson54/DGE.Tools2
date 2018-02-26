#I use this to build and locally install the package while testing.
#Not for use by end Users!

x = getwd()
setwd ("~/R/lib/pkgsrc/DGE.Tools2/")
library(devtools)
document()
load_all()
build_vignettes()
build()
pkg = "~/R/lib/pkgsrc/DGE.Tools2_0.9.36.tar.gz"
install.packages(pkg, repos=NULL, type="source")
setwd(x)

#install from Git
#After pushing to git...
install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools2")
install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools2", repos=BiocInstaller::biocinstallRepos())

#for dev
install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools", branch="dev", repos=BiocInstaller::biocinstallRepos())


#update BRAN
library(bmsPackageTools)
addPackageToRepo()
updateBmsRepo()

function todo

	extractCol.R  (rewrite to check for matching rownames/order)







