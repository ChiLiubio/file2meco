#' Transform 'microtable' object of 'microeco' package to the 'phyloseq' object of 'phyloseq' package.
#'
#' @param dataset a microtable object.
#' @return phyloseq object.
#' @examples
#' \dontrun{
#' library(microeco)
#' data("dataset")
#' meco2phyloseq(dataset)
#' }
#' @export
meco2phyloseq <- function(dataset){
	otu_table_trans <- dataset$otu_table
	sample_table_trans <- dataset$sample_table
	tax_table_trans <- dataset$tax_table
	phylo_tree_trans <- dataset$phylo_tree
	seq_trans <- dataset$rep_fasta
	if(!is.null(seq_trans)){
		if(!inherits(seq_trans, "DNAStringSet")){
			message("Skip rep_fasta conversion as it is not DNAStringSet class!")
			seq_trans <- NULL
		}
	}
	# OTU table for phyloseq
	otu_table_trans <- phyloseq::otu_table(otu_table_trans, taxa_are_rows = TRUE)
	sampledata <- phyloseq::sample_data(sample_table_trans)
	rownames(sampledata) <- rownames(sample_table_trans)
	# Taxonomic table for phyloseq
	if(!is.null(tax_table_trans)){
		tax_table_trans <- phyloseq::tax_table(as.matrix(tax_table_trans))
	}
	physeq <- phyloseq::phyloseq(otu_table_trans, tax_table_trans, sampledata, phylo_tree_trans, seq_trans)
	physeq
}


#' Transform the 'phyloseq' object of 'phyloseq' package to 'microtable' object of 'microeco' package.
#'
#' @param physeq a phyloseq object.
#' @param ... parameter passed to microtable$new function of microeco package, such as auto_tidy parameter.
#' @return microtable object.
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data("GlobalPatterns")
#' phyloseq2meco(GlobalPatterns)
#' }
#' @export
phyloseq2meco <- function(physeq, ...){
	if(physeq@otu_table@taxa_are_rows){
		otu_table_trans <- as.data.frame(physeq@otu_table@.Data, check.names = FALSE, stringsAsFactors = FALSE)
	}else{
		otu_table_trans <- as.data.frame(t(physeq@otu_table@.Data), check.names = FALSE, stringsAsFactors = FALSE)
	}
	sample_table_trans <- data.frame(phyloseq::sample_data(physeq), check.names = FALSE, stringsAsFactors = FALSE)
	tax_table_trans <- as.data.frame(physeq@tax_table@.Data, check.names = FALSE, stringsAsFactors = FALSE)
	tax_table_trans %<>% tidy_taxonomy
	phylo_tree_trans <- physeq@phy_tree
	seq_trans <- physeq@refseq

	dataset <- microtable$new(sample_table = sample_table_trans, otu_table = otu_table_trans, 
		tax_table = tax_table_trans, phylo_tree = phylo_tree_trans, rep_fasta = seq_trans, ...)
	dataset
}
