#' Transform 'QIIME' results to 'microtable' object.
#'
#' @description
#' Transform 'QIIME' results to microtable object. 
#' The QIIME results refer in particular to the files of qiime1 software.
#' 
#' @param feature_table the otu table generated from 'QIIME'. Taxonomic information should be in the end of the file.
#' @param sample_table default NULL; sample metadata table; If provided, must be one of the several types of formats: 
#'   1) comma seperated file with the suffix csv or tab seperated file with suffix tsv/txt; 
#'   2) Excel type file with the suffix xlsx or xls; require \code{readxl} package to be installed; 
#'   3) \code{data.frame} object from R.
#' @param match_table default NULL; a two column table used to replace the sample names in feature table; Must be two columns without column names;
#'   The first column must be raw sample names same with those in feature table, 
#'   the second column must be new sample names same with the rownames in sample_table; Please also see the example files.
#'   If provided, must be one of the several types of formats: \cr
#'   1) comma seperated file with the suffix csv or tab seperated file with suffix tsv/txt; \cr
#'   2) Excel type file with the suffix xlsx or xls; require \code{readxl} package to be installed; \cr
#'   3) \code{data.frame} object from R.
#' @param phylo_tree default NULL; the phylogenetic tree; generally, a file with suffix "tre".
#' @param rep_fasta default NULL; the representative sequences; a fasta file, generally with suffix "fasta" or "fna" or "fa".
#' @param ... parameter passed to microtable$new function of microeco package, such as \code{auto_tidy} parameter.
#' @return \code{microtable} object.
#' @examples
#' \dontrun{
#' # use the raw data files stored inside the package
#' otu_file_path <- system.file("extdata", "otu_table_raw.txt", package="file2meco")
#' sample_file_path <- system.file("extdata", "sample_info.csv", package="file2meco")
#' phylo_file_path <- system.file("extdata", "rep_phylo.tre", package="file2meco")
#' rep_fasta_path <- system.file("extdata", "rep.fna", package="file2meco")
#' qiime1meco(otu_file_path, sample_table = sample_file_path)
#' qiime1meco(otu_file_path, sample_table = sample_file_path, 
#'   phylo_tree = phylo_file_path)
#' qiime1meco(otu_file_path, sample_table = sample_file_path, 
#'   phylo_tree = phylo_file_path, rep_fasta = rep_fasta_path)
#' }
#' @export
qiime1meco <- function(feature_table, sample_table = NULL, match_table = NULL, phylo_tree = NULL, rep_fasta = NULL, ...){
	# check whether there is a commented line
	tryread <- readLines(feature_table)
	comlines <- sum(grepl("^#", tryread))
	commented <- ifelse(comlines > 1, TRUE, FALSE)
	
	# read and parse otu_table
	otu_raw_table <- read_qiime_otu_table(feature_table, commented = commented)
	# obtain the otu table data.frame
	feature_table <- as.data.frame(otu_raw_table[[3]])
	colnames(feature_table) <- unlist(otu_raw_table[[1]])
	# obtain the taxonomic table  data.frame
	taxonomy_table_1 <- as.data.frame(split_assignments(unlist(otu_raw_table[[4]])))
	# make the taxonomic table clean, this is very important
	taxonomy_table_1 %<>% tidy_taxonomy
	# first check the match_table
	if(!is.null(match_table)){
		feature_table <- check_match_table(match_table = match_table, abund_new = feature_table)
	}
	# read sample metadata table, data.frame, row.names = 1 set rownames
	if(!is.null(sample_table)){
		sample_table <- check_sample_table(sample_table = sample_table)
	}
	# read the phylogenetic tree
	if(!is.null(phylo_tree)){
		phylo_tree <- read.tree(phylo_tree)
	}
	if(!is.null(rep_fasta)){
		rep_fasta <- seqinr::read.fasta(rep_fasta)
	}
	# create a microtable object
	dataset <- microtable$new(sample_table = sample_table, otu_table = feature_table, tax_table = taxonomy_table_1, phylo_tree = phylo_tree, rep_fasta = rep_fasta, ...)
	dataset
}



####################################################################
# The following functions are copied from the qiimer package https://cran.r-project.org/src/contrib/Archive/qiimer/
# Because the qiimer is not in CRAN now. To install it is not easy.
# All rights reserved.



# Parse a 'QIIME' OTU table file in "calssic" format.
#
# param filepath Path to OTU table file.
# param commented TRUE if the header line is preceeded by an additional
#   comment line, otherwise FALSE.  This is usually the case for OTU
#   tables generated with 'QIIME', so we default to TRUE.
#   param metadata TRUE if the OTU table contains a metadata column, otherwise
#   FALSE.  The metadata column usually contains taxonomic assignments, and
#   must be located on the right-hand side of the table.
# return A list with four attributes: sample_ids, otu_ids, counts, and 
#   metadata, a data structure similar to that returned by the python 
#   function `qiime.parse.parse_otu_table`.  The sample_ids, otu_ids, and
#   metadata attributes are character vectors.  The counts attribute is a
#   matrix with one column per sample_id and one row per otu_id.
read_qiime_otu_table <- function(filepath, commented=TRUE, metadata=TRUE) {
  f <- file(filepath, "rt")
  header_line <- readLines(f, n=1)
  if (commented) {
    header_line <- readLines(f, n=1)
  }
  col_names <- strsplit(header_line, "\t")[[1]]

  col_classes <- rep("numeric", times=length(col_names))
  col_classes[1] <- "character"
  if (metadata) {
    col_classes[length(col_classes)] <- "character"
  }

  full_otu_table <- read.table(
    f, col.names=col_names, colClasses=col_classes, sep="\t", 
    quote="", as.is=TRUE, header=FALSE)
  close(f)

  data_cols <- if (metadata) {
    2:(length(col_names) - 1) 
  } else {
    2:length(col_names)
  } 

  sample_ids <- col_names[data_cols]
  otu_ids <- as.character(full_otu_table[,1])

  counts <- as.matrix(full_otu_table[,data_cols])
  rownames(counts) <- otu_ids

  if (metadata) {
    metadata_vals <- as.character(full_otu_table[,length(col_names)])
    names(metadata_vals) <- otu_ids
  } else {
    metadata_vals <- NULL
  }
    
  list(
    sample_ids=sample_ids, otu_ids=otu_ids, counts=counts,
    metadata=metadata_vals)
}

# Standard taxonomic ranks.
taxonomic_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Split taxonomic assignment strings
#
# param assignments Character vector of taxonomic assignments.
# param ranks Character vector of taxonomic ranks, used as column names in the
#   resultant data frame.
# param split Pattern on which to split taxa in assignment strings.
# param ... Additional parameters are passed to the \code{strsplit} function.
# return A data frame of taxonomic assignments.
split_assignments <- function(assignments, ranks=taxonomic_ranks, 
  split="; ", ...) {
  a <- strsplit(as.character(assignments), split, ...)
  max_ranks <- max(sapply(a, length))
  a <- lapply(a, function (x) {
    fill_length <- max_ranks - length(x)
    c(x, rep(NA, fill_length))
  })
  a <- as.data.frame(do.call(rbind, a))
  colnames(a) <- ranks[1:ncol(a)]
  if (!is.null(names(assignments))) {
    rownames(a) <- names(assignments)
  }
  a
}

