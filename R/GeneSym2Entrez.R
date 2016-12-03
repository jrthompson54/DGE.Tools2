### Function GeneSym2Entrez ###
#' Function GeneSym2Entrez:  Map Official Gene Symbols to Entrez GeneIDs
#'
#' Maps Gene Symbols to Entrez GeneID using the org.Xx.eg.db packages
#' for Human, Mouse and Rat. Aggregates duplicate mappings (genesym mapped to 2
#' or more Entrez IDs) to avoid row expansion.
#'
#' Note: currently aggregate mappings are lost due to exact matching.  Need
#' to implement a second pass for unmapped genesymbols to use partial matching.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID
#'
#' @param DF a dataframe with a gene symbol column
#' @param Species One of "Human", "Mouse", "Rat", "Dog" or NULL.
#' @param GeneSym  Column name of the gene symbol column to map
#'
#' @return DataFrame with Entrez GeneID column added
#'
#' @examples
#' MyDataframe = GeneSym2Entrez (MyDataframe, Species = "Human", GeneName = "GeneSym")
#'
#' @import stringi dplyr AnnotationDbi
#'
#' @export
GeneSym2Entrez = function(DF, Species=NULL, GeneSym="GeneName"){
  #JRT 20Oct2015
  #Map official gene symbols (not alias) to EntrezIDs
  #Add an EntrezID column and return the input DF with the new ID column
  #preserve rownames on the returned DF
  #aggregate duplicate mappings to avoid row expansion
  #
  #Parameters:
  #  DF = dataframe with a GeneSymbol column
  #  Species = Human|Mouse|Rat|Dog (Default = Human)
  #  GeneSym = name of the genesymbol column (Default = "GeneName" which is the
  #  default name for the GeneSymbol column in Omicsoft RefGene Annotation)

  FVersion = "05Jan2016"

  if (is.null(Species)){
    stop("Species argument is required [\"human\", \"mouse\" , \"rat\" or \"dog\"]")
  }
  if (is.null(DF[[GeneSym]])){
    stop(paste("Column", GeneSym, "not found in dataframe."))
  }

  #Get the right db
  if (tolower(Species) == "human"){
    e2s = AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
  } else if (tolower(Species) == "mouse") {
    e2s = AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egSYMBOL)
  } else if (tolower(Species) == "rat") {
    e2s = AnnotationDbi::toTable(org.Rn.eg.db::org.Rn.egSYMBOL)
  } else if (tolower(Species) == "dog") {
    e2s = AnnotationDbi::toTable(org.Cf.eg.db::org.Cf.egSYMBOL)
  }

  #aggregate duplicate mappings to avoid row expansion
  #results in one row per entrezID with multiple mapping gene symbols
  #in a comma separated list.
  e2s = aggregate(gene_id ~ symbol, data = e2s, paste, collapse=", ")
  colnames(e2s)=c(GeneSym, "EntrezID")
  #may need to upcase the GeneSym and then leftjoin on the upcased GeneSym

  #put rownames in the table
  DF$zzRownames = rownames(DF)
  rowcount = nrow(DF)
  #Left join the EntrezID on GeneSym
  DF = dplyr::left_join (DF, e2s)
  newrowcount = nrow(DF)
  #Put the rownames back
  rownames(DF) = DF$zzRownames
  #remove the extra column.
  DF$zzRownames=NULL

  if (newrowcount > rowcount) {#this should never happen given the aggregation we do above
    Warning("Warning: EntrezID addition expanded the rowcount!")
  }

  #todo
  #Add a second pass here.  Where EntrezID is still null, use partial matching of
  #the GeneSym in DF to the comma separated list of GeneSymbols for EntrezIDs
  #that are mapped to multiple genesymbols

  return(DF)
}

# Notes:
#
# org.Hs.eg.db returns all upcase symbols
#
# org.Mm.eg.db returns genesymbols with mixed case.
# Most but not all start with an uppercase letter.
#
# about 290 mouse symbols map to more than one Entrez ID.
# > e2s = AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egSYMBOL)
# > nrow(e2s)
# [1] 73347
# > gsunique = unique(e2s$symbol)
# > length(gsunique)
# [1] 73057
#
# no worse if we uppercase the symbol:
#   > gsup = toupper(e2s$symbol)
# > gsunique = unique(e2s$symbol)
# > length(gsunique)
# [1] 73057
#
# Omicsoft annotation is all uppercase genesym for human
# Omicsoft annotation is first letter upcase genesym for mouse
