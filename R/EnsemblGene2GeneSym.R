### Function  EnsemblGene2GeneSym ###
#' Function  EnsemblGene2GeneSym
#'
#' Takes a dataframe from topTable or topTreat and maps Ensembl GeneIDs to Gene
#' Symbols using the Annotables package.  Technically it should work on
#' any dataframe with Ensembl Gene IDs as rowlabels.  The added columns to the
#' include "symbol", "biotype" and "description".
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Ensembl, GeneID
#'
#' @param DF a dataframe with a gene symbol column
#' @param Build One of grch38, grch37, grcm38, rnor6.
#' @param MergeAll If TRUE, will merge all annotation columns from annotables.
#'   If FALSE (default), will only merge symbol, biotype and description columns.
#'
#' @return Dataframe with Entez GeneID, biotype and desscripton columns added
#'
#' @examples
#' #Annotate one topTable dataframe
#' MyDataframe = EnsemblGene2GeneSym (MyDataframe, Build = "grch38")
#'
#' #same with less typing
#' MyDataframe %<>% EnsemblGene2GeneSym("grch38")
#'
#' #Annotate a whole bunch of topTable dataframes in a contrast list
#' #produced by runContrasts
#' for (i in 1:length(ContrastList$TopTableList)){
#'    ContrastList$TopTableList[[i]] %<>% EnsemblGene2GeneSym("grcm38")
#' }
#'
#' @import annotables dplyr
#'
#' @export
EnsemblGene2GeneSym = function(DF, Build=NULL, MergeAll = FALSE){

  #ensgene is the ensemblgene id
  #symbol: gene symbol
  #biotype: protein, lncRNA etc
  #description:

  FVersion = "9Mar2016"
  if (!requireNamespace("annotables", quietly = TRUE)) {
    stop("package \"annotables\" is needed.  Please load it from  https://github.com/stephenturner/annotables",
         call. = FALSE)
  }
  if (is.null(Build)){
    stop("Build argument is required [\"grch38\", \"grch37\", \"grcm38\" or \"rnor6\"]")
  }

  rowcount = nrow(DF)

  #Get the right db
  assign ("MyAnnotation", get(Build))
  #need to make it unique on ensgene;  There are duplicates due to 1 to many
  #mappings to Entrez
  MyAnnotation <- MyAnnotation[!duplicated(MyAnnotation[,"ensgene"]),]

  #convert row names to a column
  DF$ensgene <- rownames(DF)

  if (!MergeAll) {
   #select the desired columns to merge
   MyAnnotation %<>% select(ensgene, symbol, biotype, description)
  }

  DF %<>% left_join(MyAnnotation, by = "ensgene")

  #put the rownames back
  rownames(DF) <- DF$ensgene
  DF$ensgene <- NULL

  if (nrow(DF) > rowcount) {
    warning("Rowcount expansion occured during annotation merge.")
  }

  return(DF)
}

