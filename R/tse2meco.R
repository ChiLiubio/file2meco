#' Transform the 'TreeSummarizedExperiment' object to 'microtable' object of 'microeco' package.
#'
#' @param tse a TreeSummarizedExperiment object <doi:10.12688/f1000research.26669.2>.
#' @param ... parameter passed to microtable$new function of microeco package, such as auto_tidy parameter.
#' @return microtable object.
#' @export
tse2meco <- function(tse, ...){
	if(!inherits(tse, "TreeSummarizedExperiment")){
		stop("Currently only support TreeSummarizedExperiment class!")
	}
	
	feature_table_trans <- as.data.frame(tse@assays@data[[1]])

	sample_table_trans <- SummarizedExperiment::colData(tse)
	if(!is.null(sample_table_trans)){
		sample_table_trans <- as.data.frame(sample_table_trans)
	}
	
	tax_table_trans <- SummarizedExperiment::rowData(tse)
	if(!is.null(tax_table_trans)){
		tax_table_trans <- as.data.frame(tax_table_trans)
	}

	if(inherits(tse, "TreeSummarizedExperiment")){
		phylo_tree_trans <- tse@rowTree$phylo
		seq_trans <- tse@referenceSeq
	}else{
		phylo_tree_trans <- NULL
		seq_trans <- NULL
	}

	meco <- microtable$new(otu_table = feature_table_trans, sample_table = sample_table_trans, 
		tax_table = tax_table_trans, phylo_tree = phylo_tree_trans, rep_fasta = seq_trans, ...)
	meco
}


#' Transform 'microtable' object of 'microeco' package to the 'TreeSummarizedExperiment' object of 'TreeSummarizedExperiment' package <doi:10.12688/f1000research.26669.2>.
#'
#' @param meco a microtable object.
#' @param ... parameter passed to \code{TreeSummarizedExperiment} function of TreeSummarizedExperiment package, e.g. \code{colTree}, \code{rowNodeLab} and \code{colNodeLab}.
#' @return TreeSummarizedExperiment object.
#' @examples
#' \dontrun{
#' library(microeco)
#' data("dataset")
#' meco2tse(dataset)
#' }
#' @export
meco2tse <- function(meco, ...){
	if(!inherits(meco, "microtable")){
		stop("The input data should be microtable object!")
	}
	if(!is.null(meco$rep_fasta)){
		if(!inherits(meco$rep_fasta, "DNAStringSet")){
			repseq <- NULL
		}else{
			repseq <- meco$rep_fasta
		}
	}else{
		repseq <- NULL
	}
	tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(assays = list(Count = meco$otu_table),
		rowData = meco$tax_table,
		colData = meco$sample_table,
		rowTree = meco$phylo_tree,
		referenceSeq = repseq,
		...)
	
	tse
}
