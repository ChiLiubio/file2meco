#' Transform metagenomic classification results of 'mpa' format to 'microtable' object.
#'
#' @description
#' Transform the classification results of mpa (MetaPhlAn) format to microtable object,
#' such as MetaPhlAn and Kraken2 results. Kraken2 results can be obtained by merge_metaphlan_tables.py from MetaPhlAn or 
#' combine_mpa.py from KrakenTools (https://ccb.jhu.edu/software/krakentools/).
#' 
#' @param abund_table 'mpa' format abundance table, see the example.
#' @param sample_data default NULL; the sample metadata table, must be tab or comma seperated file, generally, a file with suffix "tsv" or "csv"..
#' @param match_table default NULL; a two column table used to replace the sample names in 'HUMAnN abundance result; Remember just two columns with no column names;
#'    The first column must be sample names used in abund_table, the second column is the new sample names, e.g. the rownames in sample_table. See the example files.
#' @return microtable object.
#' @examples
#' \donttest{
#' # use the raw data files stored inside the package
#' abund_file_path <- system.file("extdata", "example_kraken2_merge.txt", package="file2meco")
#' sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", 
#'   package="file2meco")
#' match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
#' library(microeco)
#' library(file2meco)
#' library(magrittr)
#' mpa2meco(abund_table = abund_file_path)
#' test <- mpa2meco(abund_table = abund_file_path, sample_data = sample_file_path, 
#'   match_table = match_file_path)
#' test$tidy_dataset()
#' }
#' @export
mpa2meco <- function(abund_table, sample_data = NULL, match_table = NULL){
	# first check func_data file format.
	abund_raw <- read.delim(abund_table, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
	# extract species data
	abund_new <- abund_raw[grepl("s__", rownames(abund_raw)), ]

	# generate the taxonomic table
	raw_taxonomy <- rownames(abund_new)
	# colnames(taxonomy)[1] <- "Taxon"
	# Because of parts of taxonomy are missing, we use extracting for each taxonomy
	taxonomy <- matrix(nrow = nrow(abund_new), ncol = 7)
	colnames(taxonomy) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
	rownames(taxonomy) <- raw_taxonomy
	taxonomy[, 1] <- gsub("((d|k)__.*?)\\|.*", "\\1", raw_taxonomy, ignore.case = TRUE)
	taxonomy[, 1][!grepl("(d|k)__", taxonomy[, 1])] <- ""
	taxonomy[, 2] <- gsub(".*(p__.*?)\\|.*", "\\1", raw_taxonomy, ignore.case = TRUE)
	taxonomy[, 2][!grepl("p__", taxonomy[, 2])] <- ""
	taxonomy[, 3] <- gsub(".*(c__.*?)\\|.*", "\\1", raw_taxonomy, ignore.case = TRUE)
	taxonomy[, 3][!grepl("c__", taxonomy[, 3])] <- ""
	taxonomy[, 4] <- gsub(".*(o__.*?)\\|.*", "\\1", raw_taxonomy, ignore.case = TRUE)
	taxonomy[, 4][!grepl("o__", taxonomy[, 4])] <- ""
	taxonomy[, 5] <- gsub(".*(f__.*?)\\|.*", "\\1", raw_taxonomy, ignore.case = TRUE)
	taxonomy[, 5][!grepl("f__", taxonomy[, 5])] <- ""
	taxonomy[, 6] <- gsub(".*(g__.*?)\\|.*", "\\1", raw_taxonomy, ignore.case = TRUE)
	taxonomy[, 6][!grepl("g__", taxonomy[, 6])] <- ""
	taxonomy[, 7] <- gsub(".*(s__.*?)", "\\1", raw_taxonomy, ignore.case = TRUE)
	taxonomy[, 7][!grepl("s__", taxonomy[, 7])] <- ""

	taxonomy %<>% as.data.frame(stringsAsFactors = FALSE)
	# taxonomy <- suppressWarnings(taxonomy %>% tidyr::separate(Taxon, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = tax_sep))
	tax_table <- taxonomy

	# read sample metadata table, data.frame, row.names = 1 set rownames
	if(!is.null(sample_data)){
		if(grepl("csv", sample_data)){
			sample_data <- read.csv(sample_data, row.names = 1, stringsAsFactors = FALSE)
		}else{
			sample_data <- read.delim(sample_data, row.names = 1, stringsAsFactors = FALSE)
		}
		if(!is.null(match_table)){
			if(grepl("csv", match_table)){
				match_table <- read.csv(match_table, stringsAsFactors = FALSE, header = FALSE)
			}else{
				match_table <- read.table(match_table, stringsAsFactors = FALSE, sep = "\t")
			}
			rownames(match_table) <- match_table[, 1]
			abund_new %<>% .[, rownames(match_table), drop = FALSE]
			colnames(abund_new) <- match_table[, 2]
		}
	}

	dataset <- microtable$new(otu_table = abund_new, sample_table = sample_data, tax_table = tax_table)
	dataset
}
