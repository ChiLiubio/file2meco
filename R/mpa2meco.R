#' Transform metagenomic classification results of 'mpa' format to 'microtable' object.
#'
#' @description
#' Transform the classification results of mpa (MetaPhlAn) format to microtable object,
#' such as MetaPhlAn and Kraken2/Bracken results. Kraken2/Bracken results can be obtained by merge_metaphlan_tables.py from MetaPhlAn or 
#' combine_mpa.py from KrakenTools (https://ccb.jhu.edu/software/krakentools/).
#' The algorithm of Kraken2 determines that the abundance of a taxon is not equal to the sum of abundances of taxa in its subordinate lineage.
#' So the default tables in taxa_abund of return microtable object are extracted from the abundances of raw file. 
#' It is totally different with the return taxa_abund of cal_abund function, 
#' which sums the abundances of taxa at different taxonomic levels based on the taxonomic table and the otu_table 
#' (i.e., taxa abundance table at a specified level, e.g., 's__').
#' 
#' @param feature_table 'mpa' format abundance table, see the example.
#' @param sample_table default NULL; sample metadata table; If provided, must be one of the several types of formats: 
#'   1) comma seperated file with the suffix csv or tab seperated file with suffix tsv/txt; 
#'   2) Excel type file with the suffix xlsx or xls; require \code{readxl} package to be installed; 
#'   3) \code{data.frame} object from R.
#' @param match_table default NULL; a two column table used to replace the sample names in abundance table; Must be two columns without column names;
#'    The first column must be raw sample names same with those in feature table, 
#'    the second column must be new sample names same with the rownames in sample_table; Please also see the example files.
#' @param use_level default "s__"; the prefix parsed for the otu_table and tax_table; must be one of 'd__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__' and 's__'.
#' @param rel default FALSE; Whether convert the original abundance to relative abundance in generated taxa_abund list. 
#'    If TRUE, all the data.frame objects in taxa_abund list have relative abundance (0-1).
#' @param sel_same default 1; How to select the taxonomic information in \code{tax_table} when multiple taxonomic information have same prefix (e.g., "k__Eukaryota|k__Fungi|"). 
#'    1 represents the first one. 2 denotes the second one. 3 means selecting both.
#'    If \code{sel_same = 3}, the "|" in the taxonomic information will be replaced with ":".
#'    For the \code{taxa_abund} list, both names are remained to facilitate subsequent filtering, and the "|" will be also replaced with ":".
#' @param sep default "\\t"; The separator in provided \code{feature_table}.
#' @param skip_comments default TRUE. Whether skip the comment lines starting with "#".
#' @param ... parameter passed to microtable$new function of microeco package, such as auto_tidy parameter.
#' @return microtable object.
#' @examples
#' \donttest{
#' library(microeco)
#' library(file2meco)
#' library(magrittr)
#' # use Kraken2 file stored inside the package
#' abund_file_path <- system.file("extdata", "example_kraken2_merge.txt", package="file2meco")
#' mpa2meco(abund_file_path)
#' # add sample information table
#' sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", 
#'   package="file2meco")
#' # sample names are different between abund_file_path and sample_file_path; 
#' # use a matching table to adjust them
#' match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
#' test <- mpa2meco(abund_file_path, sample_table = sample_file_path, 
#'   match_table = match_file_path, use_level = "s__", rel = TRUE)
#' # make the taxonomy standard for the following analysis
#' test$tax_table %<>% tidy_taxonomy
#' test$tidy_dataset()
#' # calculate taxa_abund with specified level instead of raw kraken results
#' test1 <- clone(test)
#' test1$cal_abund(rel = TRUE)
#' identical(test$taxa_abund$Kingdom, test1$taxa_abund$Kingdom)
#' }
#' @export
mpa2meco <- function(feature_table, sample_table = NULL, match_table = NULL, use_level = "s__", rel = FALSE, sel_same = 1, sep = "\t", skip_comments = TRUE, ...){
	abund_raw <- readLines(feature_table)
	# check the comment lines
	if(skip_comments){
		start_line <- which(!grepl("^#", abund_raw))[1]
		abund_raw <- abund_raw[start_line:length(abund_raw)]
	}
	header_line <- unlist(strsplit(abund_raw[1], sep))
	if(length(header_line) == 1){
		stop("Please check the separator of input feature_table! Try to reset the parameter sep!")
	}
	total_abund <- sapply(abund_raw[-1], function(x){unlist(strsplit(x, sep))}) %>% 
		t %>% 
		as.data.frame %>%
		`colnames<-`(header_line) %>%
		`row.names<-`(.[, 1]) %>%
		.[, -1, drop = FALSE] %>%
		microeco::dropallfactors(unfac2num = TRUE)

	if(! use_level %in% c('d__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')){
		stop("use_level must be one of 'd__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__' and 's__'!")
	}
	
	replace_level <- c('k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')
	
	# Because parts of taxonomy are missing, we extract taxonomy for each one
	# determine the prefix of Kingdom or Domain, e.g., Bacteria, Archaea or Fungi
	if(use_level == "d__"){
		use_level = "k__"
	}
	if(any(grepl("d__", rownames(total_abund)))){
		message("Replace the prefix d__ with k__ to unify the taxonomic information ...")
		rownames(total_abund) %<>% gsub("d__", "k__", .)
	}
	all_taxonomic_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
	# extract data to create microtable object
	abund_table_taxa <- total_abund[grepl(paste0(use_level, "(?!.*\\|).*"), rownames(total_abund), perl = TRUE), , drop = FALSE]
	# generate taxonomic table
	raw_taxonomy <- rownames(abund_table_taxa)

	col_number <- which(replace_level %in% use_level)
	message("Generate otu_table at ", all_taxonomic_levels[col_number], " level ...")

	taxonomy <- matrix(nrow = nrow(abund_table_taxa), ncol = col_number)
	colnames(taxonomy) <- all_taxonomic_levels[1:col_number]
	rownames(taxonomy) <- raw_taxonomy
	for(i in 1:col_number){
		tmp <- replace_level[i]
		if(tmp == "k__"){
			if(sel_same == 1){
				taxonomy[, i] <- gsub("^(k__.*?)\\|.*", "\\1", raw_taxonomy, ignore.case = TRUE)
			}else{
				tmp_replace <- gsub(paste0(paste0(replace_level[-1], ".*"), collapse = "|"), "", raw_taxonomy, ignore.case = TRUE) %>%
					gsub("\\|$", "", .)
				if(sel_same == 2){
					taxonomy[, i] <- gsub(".*\\|", "", tmp_replace, ignore.case = TRUE)
				}else{
					if(sel_same == 3){
						taxonomy[, i] <- tmp_replace
						taxonomy[, i] %<>% gsub("\\|", ":", .)
					}else{
						stop("Parameter sel_same should be 1, 2 or 3!")
					}
				}
			}
		}else{
			taxonomy[, i] <- gsub(paste0(".*(", tmp, ".*?)(\\|.*|$)"), "\\1", raw_taxonomy, ignore.case = TRUE)
		}
		taxonomy[, i][!grepl(tmp, taxonomy[, i])] <- ""
	}
	taxonomy %<>% as.data.frame(stringsAsFactors = FALSE)
	taxonomy %<>% microeco::tidy_taxonomy()
	tax_table <- taxonomy
	message("Generate tax_table ...")
	# generate taxa_abund from raw data
	taxa_abund <- list()
	for(i in 1:col_number){
		tmp <- replace_level[i]
		tmp_table <- total_abund[grepl(paste0(tmp, "(?!.*\\|).*"), rownames(total_abund), perl = TRUE), , drop = FALSE]
		if(tmp == "k__"){
			rownames(tmp_table) %<>% gsub("\\|", ":", .)
		}
		taxa_abund[[all_taxonomic_levels[i]]] <- tmp_table
	}
	
	# first check the match_table
	if(!is.null(match_table)){
		abund_table_taxa <- check_match_table(match_table = match_table, abund_new = abund_table_taxa)
		for(i in names(taxa_abund)){
			taxa_abund[[i]] <- check_match_table(match_table = match_table, abund_new = taxa_abund[[i]])
		}
	}
	# read sample metadata table
	if(!is.null(sample_table)){
		sample_table <- check_sample_table(sample_table = sample_table)
	}
	# create microtable object
	dataset <- microtable$new(otu_table = abund_table_taxa, sample_table = sample_table, tax_table = tax_table, ...)
	message("Create the microtable object ...")
	if(rel){
		taxa_abund %<>% lapply(function(x){as.data.frame(apply(x, 2, function(y){y/sum(y)}))})
	}
	dataset$taxa_abund <- taxa_abund
	message("Generate taxa_abund list using raw taxonomic abundance of the input file ...")
	
	dataset
}

