#' Transform viromescan results to 'microtable' object.
#'
#' @description
#' Transform the results of viromescan software to microtable object.
#' The output of viromescan is single file for each sample.
#' All the results are needed to be merged and adjusted (for several chaotic taxonomy).
#' The input should be the 'count' tables at Species level, i.e. Species_level_results-Counts.txt.
#' For more details, please see the reference <DOI: 10.1186/s12864-016-2446-3>.
#' 
#' @param input_dir the input directory, containing all the result folders for each sample. Each folder should be named by the sample name.
#' @param sample_table default NULL; sample metadata table; If provided, must be one of the several types of formats:\cr
#'   1) comma seperated file with the suffix csv or tab seperated file with suffix tsv/txt;\cr
#'   2) Excel type file with the suffix xlsx or xls; require \code{readxl} package to be installed;\cr
#'   3) \code{data.frame} object from R.
#' @param match_table default NULL; a two column table used to replace the sample names in abundance table; Must be two columns without column names;
#'    The first column must be raw sample names same with those in feature table, 
#'    the second column must be new sample names same with the rownames in sample_table; Please also see the example files.
#' @param ... parameter passed to microtable$new function of microeco package, such as auto_tidy parameter.
#' @return microtable object.
#' @examples
#' \donttest{
#' library(microeco)
#' library(file2meco)
#' # use viromescan directory inside the package
#' dir_path <- system.file("extdata", "viromescan", package="file2meco")
#' d1 <- vs2meco(dir_path)
#' d1$cal_abund()
#' # d1$taxa_abund$Family is same with the percentage output of viromescan at 
#' # Family level, i.e. Family_level_results-%.txt file
#' d1$cal_abund(rel = FALSE)
#' # d1$taxa_abund$Family is same with the count output of viromescan at 
#' # Family level, i.e. Family_level_results-Counts.txt file
#' }
#' @export
vs2meco <- function(input_dir, sample_table = NULL, match_table = NULL, ...){
	# check the input
	if(!dir.exists(input_dir)){
		stop("The input_dir does not exist!")
	}
	allsamples <- list.files(input_dir)
	if(length(allsamples) == 0){
		stop("No directory found in the input_dir: ", input_dir, " !")
	}
	res_species <- NULL
	for(i in allsamples){
		# species
		path <- paste0(input_dir, "/", i, "/")
		if(!dir.exists(path)){
			stop("The sample directory: ", path, " does not exist! Please check it!")
		}
		sample_files <- list.files(path)
		if(!any(grepl("^Species_", sample_files))){
			stop("No file named Species_level_results-Counts.txt in the folder: ", path, "! Please check it!")
		}
		sample_taxon <- sample_files[grepl("^Species_", sample_files)]
		if(!any(grepl("-Counts", sample_taxon))){
			stop("No file named Species_level_results-Counts.txt in the folder: ", path, "! Please check it!")
		}
		taxon_count_path <- paste0(path, sample_taxon[grepl("-Counts", sample_taxon)])
		f1 <- read.delim(taxon_count_path, skip = 1, header = FALSE)
		colnames(f1) <- c("taxon", i)
		if(is.null(res_species)){
			res_species <- f1
		}else{
			res_species <- merge(res_species, f1, by.x = "taxon", by.y = "taxon", all = TRUE)
		}
	}
	res_species[is.na(res_species)] <- 0
	raw_taxa <- res_species[, 1]
	if(any(grepl(";$", raw_taxa))){
		message("Find some taxonomic lineage not standard! Revise them ...")
		tmp <- raw_taxa[grepl(";$", raw_taxa)]
		if(any(grepl("^Emaravirus", tmp))){
			tmp[grepl("Emaravirus", tmp)] <- paste0("Fimoviridae;", tmp[grepl("Emaravirus", tmp)])
			tmp[grepl("Emaravirus", tmp)] <- gsub(";$", "", tmp[grepl("Emaravirus", tmp)])
		}
		if(any(grepl("^Tenuivirus", tmp))){
			tmp[grepl("Tenuivirus", tmp)] <- paste0("Phenuiviridae;", tmp[grepl("Tenuivirus", tmp)])
			tmp[grepl("Tenuivirus", tmp)] <- gsub(";$", "", tmp[grepl("Tenuivirus", tmp)])
		}
		if(any(grepl(";$", tmp))){
			tmp[grepl(";$", tmp)] <- paste0("Unclassified;", tmp[grepl(";$", tmp)])
			tmp[grepl(";$", tmp)] <- gsub(";$", "", tmp[grepl(";$", tmp)])
		}
		raw_taxa[grepl(";$", raw_taxa)] <- tmp
		res_species[, 1] <- raw_taxa
	}
	rownames(res_species) <- res_species[, 1]
	res_species <- res_species[, -1, drop = FALSE]
	abund_table_taxa <- res_species

	tax_table <- sapply(rownames(abund_table_taxa), function(x){unlist(strsplit(x, ";"))}) %>% 
			t %>% 
			as.data.frame %>%
			microeco::dropallfactors(unfac2num = TRUE)
	colnames(tax_table) <- c("Family", "Genus", "Species")
	tax_table %<>% microeco::tidy_taxonomy()
	message("Generate tax_table ...")
	# first check the match_table
	if(!is.null(match_table)){
		abund_table_taxa <- check_match_table(match_table = match_table, abund_new = abund_table_taxa)
	}
	# read and check metadata table
	if(!is.null(sample_table)){
		sample_table <- check_sample_table(sample_table = sample_table)
	}
	meco_vs <- microtable$new(otu_table = abund_table_taxa, sample_table = sample_table, tax_table = tax_table, ...)
	message("Create the microtable object ...")
	meco_vs$tidy_dataset()
	meco_vs
}
