#' @import ape
#' @import microeco
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom utils data
#' @importFrom utils read.csv 
#' @importFrom utils read.delim
#' @importFrom utils read.table
#' @importFrom yaml read_yaml
#' @importFrom rhdf5 h5read
#' @importFrom Matrix sparseMatrix
#' @importFrom tidyr separate
#' @importFrom dplyr if_else
NULL

# define globalVariables
utils::globalVariables(c(".", "CHOCOPhlAn_taxonomy", "MetaCyc_pathway_map", "Tax4Fun2_KEGG", "ncyc_map", "pcyc_map", "Taxon", "unzip", "as.dist"))
