\name{GeneData}
\alias{GeneData}
\docType{data}
\title{ 
GeneData: List of List describing how to import text files
}
\description{
Defines the data files to import into a SummarizedExperiment.
Each file is described by a list with defines the fieldname, the type of data
and the text file containing the data.  The predefined GeneData item
is preconfigured for datafiles from the Omicsoft pipeline.  Modify GeneData
file values to use data from a different pipeline.
}
\usage{data("GeneData")}
\format{
  The format is:
List of 5
 $ :List of 3
  ..$ name: chr "Gene.Annotation"
  ..$ file: chr "RNA-Seq.Count.Annotation.txt"
  ..$ type: chr "rowRanges"
 $ :List of 3
  ..$ name: chr "Design"
  ..$ file: chr "RNA-Seq.Design.txt"
  ..$ type: chr "colData"
 $ :List of 3
  ..$ name: chr "Gene.Counts"
  ..$ file: chr "RNA-Seq.Count.Table.txt"
  ..$ type: chr "Assay"
 $ :List of 3
  ..$ name: chr "Gene.FPKM"
  ..$ file: chr "RNA-Seq.FPKM.Table.txt"
  ..$ type: chr "Assay"
 $ :List of 3
  ..$ name: chr "QC.Metrics"
  ..$ file: chr "RNA-Seq.QCMetrics.Table.txt"
  ..$ type: chr "metadata"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(GeneData)
## maybe str(GeneData) ; plot(GeneData) ...
}
\keyword{datasets}
