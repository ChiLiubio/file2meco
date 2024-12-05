#' The CHOCOPhlAn_taxonomy data
#'
#' The CHOCOPhlAn_taxonomy data is used for the parsing the 'HUMAaN' metagenomic results and add the taxonomy hierarchical information to the 'tax_table'.
#'
#' @docType data
#' @keywords data.frame
#' @name CHOCOPhlAn_taxonomy
#' @usage data(CHOCOPhlAn_taxonomy)
NULL

#' The ncyc_map data
#'
#' The ncyc_map data is used for the parsing the 'NCycDB' metagenomic results and add the N cycle pathway information to the 'tax_table' of 'microtable' object.
#'
#' @docType data
#' @keywords data.frame
#' @name ncyc_map
#' @usage data(ncyc_map)
NULL

#' The pcyc_map data
#'
#' The pcyc_map data is used for the parsing the 'PCycDB' metagenomic results and add the P cycle pathway information to the 'tax_table' of 'microtable' object.
#'
#' @docType data
#' @keywords data.frame
#' @name pcyc_map
#' @usage data(pcyc_map)
NULL

#' The MetaCyc_pathway_map data
#'
#' The MetaCyc_pathway_map data is a manually curated hierarchical structure data of 'MetaCyc' pathways. 
#' It is used for parsing the 'HUMAaN' metagenomic abundance table associated with 'MetaCyc' database.
#' Currently, only superclass 1, 2 and the pathway are used in this data.
#' Some metabolic pathways may have multiple classification informations at specific levels, 
#' which are connected here with the symbol "&&". 
#' Therefore, when filtering in this table with the name of a pathway, please use regular expressions instead of direct filtering operations.
#'
#' @docType data
#' @keywords data.frame
#' @name MetaCyc_pathway_map
#' @usage data(MetaCyc_pathway_map)
NULL
