#' Transform 'Ncyc' metagenomic abundance to 'microtable' object.
#'
#' @description
#' Transform 'Ncyc' metagenomic abundance to microtable object. Reference: Qichao et al. (2019) <doi: 10.1093/bioinformatics/bty741>.
#' @param abund_table 'Ncyc' software output abundance table, see the example file.
#' @param sample_data default NULL; the sample metadata table; data.frame or character for the path; 
#' A file path must be tab or comma seperated file, generally, a file with suffix "tsv" or "csv".
#' @param match_table default NULL; data.frame or character for the path; should be two column table used to replace the sample names in abundance result; 
#' Remember just two columns with no column names;
#' The first column must be sample names used in abund_table, the second column is the new sample names, e.g. the rownames in sample_table. See the example files;
#' A file path must be tab or comma seperated file, e.g. a file with suffix "tsv" or "csv".
#' @param ... parameter passed to microtable$new function of microeco package, such as auto_tidy parameter.
#' @return microtable object.
#' @examples
#' \donttest{
#' # use the raw data files stored inside the package
#' abund_file_path <- system.file("extdata", "example_Ncyc_table.tsv", package="file2meco")
#' sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", 
#'   package="file2meco")
#' match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
#' library(microeco)
#' library(file2meco)
#' library(magrittr)
#' ncyc2meco(abund_table = abund_file_path)
#' test <- ncyc2meco(abund_table = abund_file_path, sample_data = sample_file_path, 
#'   match_table = match_file_path)
#' test$tidy_dataset()
#' # use split_group = TRUE to calculate the pathway abundance with multipe map correspondance
#' test$cal_abund(select_cols = 1:2, rel = TRUE, split_group = TRUE, split_column = "Pathway")
#' test$taxa_abund$Pathway %<>% .[!grepl("unclass", rownames(.)), ]
#' test1 <- trans_abund$new(test, taxrank = "Pathway")
#' test1$plot_bar(bar_type = "notfull")
#' # for gene abundance, no splitting on the Pathway
#' test$cal_abund(select_cols = 1:2, rel = TRUE, split_group = FALSE)
#' test$taxa_abund$Gene %<>% .[!grepl("unclass", rownames(.)), ]
#' test1 <- trans_abund$new(test, taxrank = "Gene")
#' test1$plot_bar(bar_type = "notfull")
#' }
#' @export
ncyc2meco <- function(abund_table, sample_data = NULL, match_table = NULL, ...){
	# first check func_data file format.
	abund_raw <- read.delim(abund_table, row.names = 1, comment.char = "#", stringsAsFactors = FALSE)

	seq_num <- readLines(abund_table)[1]
	seq_num <- as.numeric(gsub(".*\\s(\\d+)$",  "\\1", seq_num))
	
	# recalculate the abundance for unclassified
	unclassified <- seq_num - apply(abund_raw, 2, sum)
	abund_new <- rbind.data.frame(abund_raw, unclassified = unclassified)

	# first check the match_table
	if(!is.null(match_table)){
		abund_new <- check_match_table(match_table = match_table, abund_new = abund_new)
	}
	# read sample metadata table
	if(!is.null(sample_data)){
		sample_data <- check_sample_table(sample_data = sample_data)
	}
	
	data("ncyc_map", envir=environment())
	dataset <- microtable$new(otu_table = abund_new, tax_table = ncyc_map, sample_table = sample_data, ...)
	dataset
}
