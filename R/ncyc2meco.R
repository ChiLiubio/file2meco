#' Transform 'Ncyc' metagenomic abundance to 'microtable' object.
#'
#' @description
#' Transform 'Ncyc' metagenomic abundance to microtable object. Reference: Qichao et al. (2019) <doi: 10.1093/bioinformatics/bty741>.
#' @param abund_table 'Ncyc' software output abundance table, see the example file.
#' @param sample_data default NULL; the sample metadata table, must be tab or comma seperated file, generally, a file with suffix "tsv" or "csv"..
#' @param match_table default NULL; a two column table used to replace the sample names in HUMAnN abundance result; Remember just two columns with no column names;
#'    The first column must be sample names used in abund_table, the second column is the new sample names, e.g. the rownames in sample_table. See the example files.
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
ncyc2meco <- function(abund_table, sample_data = NULL, match_table = NULL){
	# first check func_data file format.
	abund_raw <- read.delim(abund_table, row.names = 1, comment.char = "#", stringsAsFactors = FALSE)

	seq_num <- readLines(abund_table)[1]
	seq_num <- as.numeric(gsub(".*\\s(\\d+)$",  "\\1", seq_num))
	
	# recalculate the abundance for unclassified
	unclassified <- seq_num - apply(abund_raw, 2, sum)
	abund_new <- rbind.data.frame(abund_raw, unclassified = unclassified)
	
	# read sample metadata table, data.frame, row.names = 1 set rownames
	if(!is.null(sample_data)){
		if(grepl("csv", sample_data)){
			sample_data <- read.csv(sample_data, row.names = 1, stringsAsFactors = FALSE)
		}else{
			sample_data <- read.delim(sample_data, row.names = 1, stringsAsFactors = FALSE)
		}
		if(!is.null(match_table)){
			if(grepl("csv", match_table)){
				match_table <- read.csv(sample_data, stringsAsFactors = FALSE, header = FALSE)
			}else{
				match_table <- read.table(match_table, stringsAsFactors = FALSE, sep = "\t")
			}
			rownames(match_table) <- match_table[, 1]
			abund_new %<>% .[, rownames(match_table), drop = FALSE]
			colnames(abund_new) <- match_table[, 2]
		}
	}
	data("ncyc_map", envir=environment())
	dataset <- microtable$new(otu_table = abund_new, tax_table = ncyc_map, sample_table = sample_data)
	dataset
}
