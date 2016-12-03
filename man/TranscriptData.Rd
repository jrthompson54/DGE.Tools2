\name{TranscriptData}
\alias{TranscriptData}
\docType{data}
\title{
TranscriptData: List of List describing how to import text files
}
\description{
Defines the data files to import into a SummarizedExperiment.
Each file is described by a list with defines the fieldname, the type of data
and the text file containing the data.  The predefined TranscriptData item
is preconfigured for datafiles from the Omicsoft pipeline.  Modify TranscriptData
file values to use data from a different pipeline.
}
\usage{data("TranscriptData")}
\format{
  The format is:
List of 5
 $ :List of 3
  ..$ name: chr "Transcript.Annotation"
  ..$ file: chr "RNA-Seq.Transcript_Count.Annotation.txt"
  ..$ type: chr "rowRanges"
 $ :List of 3
  ..$ name: chr "Design"
  ..$ file: chr "RNA-Seq.Design.txt"
  ..$ type: chr "colData"
 $ :List of 3
  ..$ name: chr "Transcript.Counts"
  ..$ file: chr "RNA-Seq.Transcript_Count.Table.txt"
  ..$ type: chr "Assay"
 $ :List of 3
  ..$ name: chr "Transcript.FPKM"
  ..$ file: chr "RNA-Seq.Transcript_FPKM.Table.txt"
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
data(TranscriptData)
## maybe str(TranscriptData) ; plot(TranscriptData) ...
}
\keyword{datasets}
