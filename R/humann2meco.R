#' Transform 'HUMAnN' metagenomic results to 'microtable' object.
#'
#' @description
#' Transform 'HUMAnN' metagenomic results to microtable object, reference: Franzosa et al. (2018) <doi:10.1038/s41592-018-0176-y>.
#' @param feature_table file path of 'HUMAnN' output abundance table; Please see the example.
#' @param db default "MetaCyc"; one of "MetaCyc", "KEGG" or "gene"; "MetaCyc" or "KEGG" means the input feature table is pathway abundance.
#'   "gene" represents the abundance of genes, such as 'eggNOG', 'KO' and 'EC'. 
#'   When using "gene", the generated tax_table has only taxonomic lineages and gene name, no higher functional levels.
#' @param sample_table default NULL; sample metadata table; If provided, must be one of the several types of formats: \cr
#'   1) comma seperated file with the suffix csv or tab seperated file with suffix tsv or txt; \cr
#'   2) Excel type file with the suffix xlsx or xls; require \code{readxl} package to be installed; \cr
#'   3) \code{data.frame} object from R.
#' @param match_table default NULL; a two column table used to replace the sample names in feature table; Must be two columns without column names;
#'   The first column must be raw sample names same with those in feature table, 
#'   the second column must be new sample names same with the rownames in sample_table; Please also see the example files.
#'   If provided, must be one of the several types of formats: \cr
#'   1) comma seperated file with the suffix csv or tab seperated file with suffix tsv/txt; \cr
#'   2) Excel type file with the suffix xlsx or xls; require \code{readxl} package to be installed; \cr
#'   3) \code{data.frame} object from R.
#' @param ... parameter passed to \code{microtable$new} function of \code{microeco} package, such as \code{auto_tidy} parameter.
#' @return \code{microtable} object.
#' @examples
#' \donttest{
#' library(file2meco)
#' library(microeco)
#' library(magrittr)
#' sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", 
#'   package="file2meco")
#' match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
#' # MetaCyc pathway examples
#' # use the raw data files stored inside the package for MetaCyc pathway database based analysis
#' abund_file_path <- system.file("extdata", "example_HUMAnN_MetaCyc_abund.tsv", package="file2meco")
#' # the default db is "MetaCyc"
#' humann2meco(abund_file_path, db = "MetaCyc")
#' humann2meco(abund_file_path, db = "MetaCyc", sample_table = sample_file_path, 
#'   match_table = match_file_path)
#' test <- humann2meco(abund_file_path, db = "MetaCyc", sample_table = sample_file_path, 
#'   match_table = match_file_path)
#' test$tidy_dataset()
#' # rel = FALSE sum original abundance instead of relative abundance
#' test$cal_abund(select_cols = 1:3, rel = FALSE)
#' test$taxa_abund$Superclass1 %<>% .[!grepl("unclass", rownames(.)), ]
#' # use_percentage = FALSE disable percentage for relative abundance
#' test1 <- trans_abund$new(test, taxrank = "Superclass1", ntaxa = 10, use_percentage = FALSE)
#' # reassign ylab title instead of default 'Relative Abundance'
#' test1$ylabname <- "Abundance (RPK)"
#' # bar_full = FALSE show original abundance instead of normalized 0-1
#' test1$plot_bar(facet = "Group", bar_full = FALSE)
#' # select both function and taxa
#' test$cal_abund(select_cols = c("Superclass1", "Phylum", "Genus"), rel = TRUE)
#' test1 <- trans_abund$new(test, taxrank = "Phylum", ntaxa = 10, delete_taxonomy_lineage = TRUE)
#' test1$plot_bar(facet = "Group")
#' test$taxa_abund$Phylum %<>% .[!grepl("unclass", rownames(.)), ]
#' test1 <- trans_abund$new(test, taxrank = "Phylum", ntaxa = 10, delete_taxonomy_lineage = FALSE)
#' test1$plot_bar(facet = "Group")
#' # functional biomarker
#' test$cal_abund(select_cols = 1:3, rel = TRUE)
#' test$taxa_abund$Superclass1 %<>% .[!grepl("unclass", rownames(.)), ]
#' test$taxa_abund$Superclass2 %<>% .[!grepl("unclass", rownames(.)), ]
#' test$taxa_abund$pathway %<>% .[!grepl("unclass", rownames(.)), ]
#' test1 <- trans_diff$new(test, method = "lefse", group = "Group")
#' test1$plot_diff_bar(use_number = 1:20)
#' # taxa biomarker
#' test$cal_abund(select_cols = 4:9, rel = TRUE)
#' test$taxa_abund$Phylum %<>% .[!grepl("unclass", rownames(.)), ]
#' test1 <- trans_diff$new(test, method = "lefse", group = "Group", p_adjust_method = "none")
#' test1$plot_diff_bar(threshold = 2)
#' #############################################################
#' # KEGG pathway examples
#' abund_file_path <- system.file("extdata", "example_HUMAnN_KEGG_abund.tsv", package="file2meco")
#' humann2meco(abund_file_path, db = "KEGG")
#' test <- humann2meco(abund_file_path, db = "KEGG", 
#'   sample_table = sample_file_path, match_table = match_file_path)
#' test$tax_table %<>% subset(Level.1 != "unclassified")
#' test$tidy_dataset()
#' test$cal_abund(select_cols = 1:3, rel = FALSE)
#' # use_percentage = FALSE disable percentage for relative abundance
#' test1 <- trans_abund$new(test, taxrank = "Level.2", ntaxa = 10, use_percentage = FALSE)
#' # or use ggplot2::ylab to change ylab title
#' test1$ylabname <- "Abundance (RPK)"
#' test1$plot_bar(facet = "Group", bar_full = FALSE)
#' # select both function and taxa
#' test$cal_abund(select_cols = c("Level.1", "Phylum", "Genus"), rel = TRUE)
#' test1 <- trans_abund$new(test, taxrank = "Phylum", ntaxa = 10, delete_taxonomy_lineage = FALSE)
#' test1$plot_bar(facet = "Group")
#' # functional biomarker
#' test$cal_abund(select_cols = 1:3, rel = TRUE)
#' test1 <- trans_diff$new(test, method = "lefse", group = "Group")
#' test1$plot_diff_bar(threshold = 3)
#' # taxa biomarker
#' test$cal_abund(select_cols = 4:9, rel = TRUE)
#' test1 <- trans_diff$new(test, method = "lefse", group = "Group", p_adjust_method = "none")
#' test1$plot_diff_bar(threshold = 2)
#' }
#' @export
humann2meco <- function(feature_table, db = c("MetaCyc", "KEGG", "gene")[1], sample_table = NULL, match_table = NULL, ...){
	db <- match.arg(db, c("MetaCyc", "KEGG", "gene"))
	# first check func_data file format.
	abund_raw <- read.delim(feature_table, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)

	# recalculate the abundance for unclassified
	abund_rawname <- rownames(abund_raw)

	abund_new <- abund_raw
	abund_rawname_func <- gsub("\\|.*$", "", abund_rawname)
	# find rows that need removed
	remove_number <- c()
	for(i in unique(abund_rawname_func)){
		filter_num <- which(abund_rawname_func == i)
		if(sum(abund_rawname_func == i) == 1){
			abund_rawname[filter_num] %<>% paste0(., "|unclassified")
		}else{
			if(sum(abund_rawname_func == i) == 2){
				# no taxa mapped, all unclassified!
				abund_new[filter_num[2], ] <- abund_new[filter_num[1], ]
			}else{
				abund_new[filter_num[length(filter_num)], ] <- abund_new[filter_num[1], ] - 
					apply(abund_new[filter_num[2:(length(filter_num) - 1)], , drop = FALSE], 2, sum)
			}
			remove_number <- c(remove_number, filter_num[1])
		}
	}
	
	# a new abund table with clear name
	rownames(abund_new) <- abund_rawname
	abund_new %<>% .[-remove_number, , drop = FALSE]
	abund_newname <- rownames(abund_new)

	data("CHOCOPhlAn_taxonomy", envir=environment())

	tax_table <- abund_newname %>%
		as.data.frame(stringsAsFactors = FALSE) %>%
		`colnames<-`(c("raw")) %>%
		tidyr::separate(col = "raw", into = c("Func", "Tax"), sep = "\\|") %>%
		tidyr::separate(col = "Tax", into = c("Genus", "Species"), sep = "\\.", fill = "right")
		
	# get the function and taxonomy lineage
	if(db == "MetaCyc"){
		data("MetaCyc_pathway_map", envir=environment())
		MetaCyc_pathway_map_use <- cbind.data.frame(rowname = rownames(MetaCyc_pathway_map), MetaCyc_pathway_map, stringsAsFactors = FALSE)
		# delete the descriptions behind the pathway names.
		tax_table$Func %<>% gsub(":.*$", "", .)
		
		tax_table <- tax_table %>%
			dplyr::left_join(., MetaCyc_pathway_map_use, by = c("Func" = "rowname")) %>%
			dplyr::left_join(., CHOCOPhlAn_taxonomy, by = c("Genus" = "Genus")) %>%
			.[, c("Superclass1", "Superclass2", "pathway", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
	}
	if(db == "KEGG"){
		data("Tax4Fun2_KEGG", envir=environment(), package = "microeco")
		ko_mapping_file <- Tax4Fun2_KEGG$ptw_desc
		ko_mapping_file$pathway <- rownames(ko_mapping_file)

		tax_table <- tax_table %>%
			dplyr::left_join(., ko_mapping_file, by = c("Func" = "pathway")) %>%
			dplyr::left_join(., CHOCOPhlAn_taxonomy, by = c("Genus" = "Genus")) %>%
			.[, c("Level.1", "Level.2", "Level.3", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
	}
	if(db == "gene"){
		tax_table <- tax_table %>%
			dplyr::left_join(., CHOCOPhlAn_taxonomy, by = c("Genus" = "Genus")) %>%
			.[, c("Func", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
	}
	
	# the CHOCOPhlAn_taxonomy data has been checked manually!
	if(nrow(tax_table) != length(abund_newname)){
		stop("Bad taxonomic lineage exists in the dataset! Please contact the maintainer!")
	}else{
		rownames(tax_table) <- abund_newname
	}
	tax_table[is.na(tax_table)] <- "unclassified"
	
	# first check the match_table
	if(!is.null(match_table)){
		abund_new <- check_match_table(match_table = match_table, abund_new = abund_new)
	}
	# read sample metadata table
	if(!is.null(sample_table)){
		sample_table <- check_sample_table(sample_table = sample_table)
	}
	
	dataset <- microtable$new(otu_table = abund_new, sample_table = sample_table, tax_table = tax_table, ...)
	dataset
}
