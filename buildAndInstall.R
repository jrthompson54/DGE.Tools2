#I use this to build and locally install the package while testing.
#Not for use by end Users!

x = getwd()
setwd ("~/R/lib/pkgsrc/DGE.Tools/")
pkg = "~/R/lib/pkgsrc/DGE.Tools_1.0.7.5.tar.gz"
library(devtools);document();load_all;build()
install.packages(pkg, repos=NULL, type="source")
setwd(x)


#install from Git
#After pushing to git...
install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools")
install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools", repos=BiocInstaller::biocinstallRepos())

#for dev
install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools", branch="dev", repos=BiocInstaller::biocinstallRepos())

