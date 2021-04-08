#' Transform HUMAnN metagenomic results to microtable object.
#'
#' @description
#' Transform HUMAnN metagenomic results to microtable object, reference: Franzosa et al. (2018) <doi:10.1038/s41592-018-0176-y>.
#' @param abund_table HUMAnN output abundance table, see the example.
#' @param db default "MetaCyc"; either "MetaCyc" or "KEGG"; the pathway database used in the abund_table file generation.
#' @param sample_data default NULL; the sample metadata table, must be tab or comma seperated file, generally, a file with suffix "tsv" or "csv"..
#' @param match_table default NULL; a two column table used to replace the sample names in HUMAnN abundance result; Remember just two columns with no column names;
#'    The first column must be sample names used in abund_table, the second column is the new sample names, e.g. the rownames in sample_table. See the example files.
#' @return microtable object.
#' @examples
#' \donttest{
#' library(file2meco)
#' library(microeco)
#' library(magrittr)
#' sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", package="file2meco")
#' match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
#' # MetaCyc pathway examples
#' # use the raw data files stored inside the package for MetaCyc pathway database based analysis
#' abund_file_path <- system.file("extdata", "example_HUMAnN_MetaCyc_abund.tsv", package="file2meco")
#' # the default db is "MetaCyc"
#' humann2meco(abund_table = abund_file_path, db = "MetaCyc")
#' humann2meco(abund_table = abund_file_path, db = "MetaCyc", sample_data = sample_file_path, match_table = match_file_path)
#' test <- humann2meco(abund_table = abund_file_path, db = "MetaCyc", sample_data = sample_file_path, match_table = match_file_path)
#' test$tidy_dataset()
#' # rel = FALSE donot use relative abundance
#' test$cal_abund(select_cols = 1:3, rel = FALSE)
#' test$taxa_abund$Superclass1 %<>% .[!grepl("unclass", rownames(.)), ]
#' test1 <- trans_abund$new(test, taxrank = "Superclass1", ntaxa = 10)
#' test1$plot_bar(facet = "Group", ylab_title = "Abundance (RPK)")
#' # select both function and taxa
#' test$cal_abund(select_cols = c("Superclass1", "Phylum", "Genus"), rel = TRUE)
#' test1 <- trans_abund$new(test, taxrank = "Phylum", ntaxa = 10, delete_part_prefix = T)
#' test1$plot_bar(facet = "Group")
#' test$taxa_abund$Phylum %<>% .[!grepl("unclass", rownames(.)), ]
#' test1 <- trans_abund$new(test, taxrank = "Phylum", ntaxa = 10, delete_part_prefix = T)
#' test1$plot_bar(facet = "Group")
#' # functional biomarker
#' test$cal_abund(select_cols = 1:3, rel = TRUE)
#' test$taxa_abund$Superclass1 %<>% .[!grepl("unclass", rownames(.)), ]
#' test$taxa_abund$Superclass2 %<>% .[!grepl("unclass", rownames(.)), ]
#' test$taxa_abund$pathway %<>% .[!grepl("unclass", rownames(.)), ]
#' test1 <- trans_diff$new(test, method = "lefse", group = "Group")
#' test1$plot_lefse_bar(use_number = 1:20)
#' # taxa biomarker
#' test$cal_abund(select_cols = 4:9, rel = TRUE)
#' test$taxa_abund$Phylum %<>% .[!grepl("unclass", rownames(.)), ]
#' test1 <- trans_diff$new(test, method = "lefse", group = "Group")
#' test1$plot_lefse_bar(LDA_score = 2)
#' #############################################################
#' # KEGG pathway examples
#' abund_file_path <- system.file("extdata", "example_HUMAnN_KEGG_abund.tsv", package="file2meco")
#' humann2meco(abund_table = abund_file_path, db = "KEGG")
#' test <- humann2meco(abund_table = abund_file_path, db = "KEGG", sample_data = sample_file_path, match_table = match_file_path)
#' test$tax_table %<>% subset(level1 != "unclassified")
#' test$tidy_dataset()
#' # rel = FALSE donot use relative abundance
#' test$cal_abund(select_cols = 1:3, rel = FALSE)
#' test1 <- trans_abund$new(test, taxrank = "level2", ntaxa = 10)
#' test1$plot_bar(facet = "Group", ylab_title = "Abundance (RPK)")
#' # select both function and taxa
#' test$cal_abund(select_cols = c("level1", "Phylum", "Genus"), rel = TRUE)
#' test1 <- trans_abund$new(test, taxrank = "Phylum", ntaxa = 10, delete_part_prefix = T)
#' test1$plot_bar(facet = "Group")
#' # functional biomarker
#' test$cal_abund(select_cols = 1:3, rel = TRUE)
#' test1 <- trans_diff$new(test, method = "lefse", group = "Group")
#' test1$plot_lefse_bar(LDA_score = 3)
#' # taxa biomarker
#' test$cal_abund(select_cols = 4:9, rel = TRUE)
#' test1 <- trans_diff$new(test, method = "lefse", group = "Group")
#' test1$plot_lefse_bar(LDA_score = 2)
#' }
#' @export
humann2meco <- function(abund_table, db = c("MetaCyc", "KEGG")[1], sample_data = NULL, match_table = NULL){
	# first check func_data file format.
	abund_raw <- read.delim(abund_table, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)

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
				abund_new[filter_num[length(filter_num)], ] <- abund_new[filter_num[1], ] - apply(abund_new[filter_num[2:(length(filter_num) - 1)], ], 2, sum)
			}
			remove_number <- c(remove_number, filter_num[1])
		}
	}
	
	# a new abund table with clear name
	rownames(abund_new) <- abund_rawname
	abund_new %<>% .[-remove_number, ]
	abund_newname <- rownames(abund_new)

	data("CHOCOPhlAn_taxonomy", envir=environment())

	tax_table <- abund_newname %>%
		as.data.frame(stringsAsFactors = FALSE) %>%
		`colnames<-`(c("raw")) %>%
		tidyr::separate(col = "raw", into = c("Func", "Tax"), sep = "\\|") %>%
		tidyr::separate(col = "Tax", into = c("Genus", "Species"), sep = "\\.", fill = "right")
		
	# get the function and taxonomy leneage
	if(grepl("MetaCyc", db, ignore.case = TRUE)){
		data("MetaCyc_pathway_map", envir=environment())
		MetaCyc_pathway_map_use <- cbind.data.frame(rowname = rownames(MetaCyc_pathway_map), MetaCyc_pathway_map, stringsAsFactors = FALSE)
		# delete the descriptions behind the pathway names.
		tax_table$Func %<>% gsub(":.*$", "", .)
		
		tax_table <- tax_table %>%
			dplyr::left_join(., MetaCyc_pathway_map_use, by = c("Func" = "rowname")) %>%
			dplyr::left_join(., CHOCOPhlAn_taxonomy, by = c("Genus" = "Genus")) %>%
			.[, c("Superclass1", "Superclass2", "pathway", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
		
	}else{
		if(grepl("KEGG", db, ignore.case = TRUE)){
			data("Tax4Fun2_KEGG", envir=environment(), package = "microeco")
			ko_mapping_file <- Tax4Fun2_KEGG$ptw_desc
			colnames(ko_mapping_file) <- c("pathway", "level3", "level2", "level1")
			ko_mapping_file%<>% .[, c(1, 4, 3, 2)]

			tax_table <- tax_table %>%
				dplyr::left_join(., ko_mapping_file, by = c("Func" = "pathway")) %>%
				dplyr::left_join(., CHOCOPhlAn_taxonomy, by = c("Genus" = "Genus")) %>%
				.[, c("level1", "level2", "level3", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

		}else{
			stop("Unknown db provided! Must be either MetaCyc or KEGG !")
		}
	}

	# the CHOCOPhlAn_taxonomy data has been checked manually! If duplicates in results, need carefully checking!
	if(nrow(tax_table) != length(abund_newname)){
		stop("Bad taxonomic lineage exists in the dataset! Please contact the maintainer!")
	}else{
		rownames(tax_table) <- abund_newname
	}
	tax_table[is.na(tax_table)] <- "unclassified"
	
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

	dataset <- microtable$new(otu_table = abund_new, sample_table = sample_data, tax_table = tax_table)
	dataset
}
