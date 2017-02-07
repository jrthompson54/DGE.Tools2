#I use this to build and locally install the package while testing.
#Not for use by end Users!

x = getwd()
setwd ("~/R/lib/pkgsrc/DGE.Tools2/")
library(devtools)
document()
load_all()
#build_vignettes()
build()
pkg = "~/R/lib/pkgsrc/DGE.Tools2_2.0.2.tar.gz"
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
	ThemePack.R
	volcanoPlot.R done
	profilePlot.R done
	comparePlot.R done
	plotPValHist.R done
	cdfPlot.R done
	obsPlot.R
	plotPvalHist

noChange needed

	OmicsoftRefGeneID2GeneSym.R
	EntrezGene2Ensembl.R
	EnsemblGene2GeneSym.R
	EnsemblGene2Entrez.R
	DGE_Unit_Conversion.R
	summarizeSigCounts.R
	GeneSym2Entrez.R
	Entrez2GeneSym.R
	JRT_heatmap.R (rewrite to use a list of parameters)
	extractCol.R  (rewrite to check for matching rownames/order)

INC1585521  Avaya softephone ticket


ThemePack, volcanoPlot, profilePlot, comparePlot, plotPValHist, cdfPlot, obsPlot.R, plotPvalHist, ggplotMDS

