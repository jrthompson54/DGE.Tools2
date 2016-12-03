### Function  EntrezGene2Ensembl ###
#' Function  EntrezGene2Ensembl
#'
#' Map Entrez GeneIDs to Ensembl GeneIDs using org.Xx.eg.db packages
#' for Human, Mouse, Rat or dog
#' 
#' One to Many handling: If the EntrezID maps to multiple EnsemblIDs, 
#' the EnsemblID field will contain a comma-separated
#' list of matching EnsemblIDs.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID
#'
#' @param DF a dataframe with a Ensembl gene ID column
#' @param Species One of Human, Mouse, Rat or Dog.
#' @param EntrezGene Column name of the Entrez GeneID column to map
#'
#' @return Dataframe with Entez GeneID column added
#'
#' @examples
#' MyDataframe = EnsemblGene2Entrez (MyDataframe, Species = "Human", EnsemblGene = "GeneID")
#'
#' @export
EntrezGene2Ensembl = function(DF, Species=NULL, EntrezGene="GeneID"){
  #  #JRT 10Aug2016 modified from Ensembl2Entrez
  #Map Ensembl GeneID to EntrezIDs
  #Add an EntrezID column and return the input DF with the new ID column
  #preserve rownames on the returned DF
  #warn if rowcount was expanded.
  #
  #Parameters:
  #  DF = dataframe with a GeneSymbol column
  #  Species = Human|Mouse|Rat (Default = Human)
  #  EnsembleGene = column name of the EnsembleGene GeneID column in your
  #  Gene/Transcript annotation data. (Default = "GeneID" which is the
  #  default name for the Ensembl GeneID column in Omicsoft Ensembl Annotation)
  #
  FVersion = "10Aug2016"

  if (is.null(Species)){
    stop("Species argument is required [\"human\", \"mouse\", \"rat\" or \"dog\"]")
  }

  #Get the right db
  if (tolower(Species) == "human"){
    #library(org.Hs.eg.db)
    et2en = AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL)
  } else if (tolower(Species) == "mouse") {
    #library(org.Mm.eg.db)
    et2en = AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egENSEMBL)
  } else if (tolower(Species) == "rat") {
    #library(org.Rn.eg.db)
    et2en = AnnotationDbi::toTable(org.Rn.eg.db::org.Rn.egENSEMBL)
  } else if (tolower(Species) == "dog") {
    #library(org.Rn.eg.db)
    et2en = AnnotationDbi::toTable(org.Cf.eg.db::org.Cf.egENSEMBL)
  }

  #aggregate multiple ensembl ids mapping to each entrezid to avoid row expansion
  #Aggregation concatenates multiple Ensembl IDs mapping to the same EntrezID (Comma sep).
  #The effectively prevents these genes from mapping through the Ensembl geneID but
  #is necessary to keep the EntrezID unique and avoid row expansion.
  #There were less than 2 dozen human genes but ~900 mouise genes that are affected by this.
  #Check the rowcount before/after aggregation to see how many.
  #Search the aggregated table for commas #in the EnsemblID field to see which
  #genes are affected by this.  The user will have to decide which of multiple 
  #mappings they wish to keep.

  #Consider using a partial matching strategy  i.e. Ensembl ID from Omicsoft
  #that's a partial match to the org.Xx.eg.db aggregated table.
  #For speed purposes, use left_join to join all the 1:1 mappings.  Then find the
  #rows that didn't map (NA in the EntrezID in the result Table) and map them with
  #partial matching strategy.

  if (is.null(DF[[EntrezGene]])) {
    stop (paste("EntrezGene field named", EntrezGene, "does not exist in the dataframe.", sep=" "))
  }

#  et2en = aggregate(gene_id ~ ensembl_id, data = et2en, paste, collapse=", ")
#  colnames(et2en)=c(EnsemblGene, "EntrezID")
  
  et2en = aggregate(ensembl_id ~ gene_id, data = et2en, paste, collapse=", ")
  colnames(et2en)=c(EntrezGene, "EnsemblGene")

  #put rownames (EntrezID) in the table
  DF$zzRownames = rownames(DF)
  rowcount = nrow(DF)
  #Left join the EntrezID
  DF = dplyr::left_join (DF, et2en)
  newrowcount = nrow(DF)
  #Put the rownames back
  rownames(DF) = DF$zzRownames
  #remove the extra column.
  DF$zzRownames=NULL
  if (newrowcount > rowcount) {
    Warning("Warning: EntrezID addition expanded the rowcount!")
  }
  return(DF)
}
