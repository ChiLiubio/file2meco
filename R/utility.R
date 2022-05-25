#' Replace the names use match table
#'
#' @param match_table default NULL; character or data.frame; matching table used.
#' @param abund_new default NULL; data.frame; the abundance table used.
#' @return new abundance table.
check_match_table <- function(match_table = NULL, abund_new = NULL){
	# read according to the input class
	if(inherits(match_table, "character")){
		if(grepl("csv", match_table)){
			match_table <- read.csv(match_table, stringsAsFactors = FALSE, header = FALSE)
		}else{
			match_table <- read.table(match_table, stringsAsFactors = FALSE, sep = "\t")
		}
	}else{
		if(! inherits(match_table, "data.frame")){
			stop("The input match_table is not data.frame class!")
		}
	}
	rownames(match_table) <- match_table[, 1]
	abund_new %<>% .[, rownames(match_table), drop = FALSE]
	colnames(abund_new) <- match_table[, 2]
	# output new abundance table
	abund_new
}

#' Read sample table
#'
#' @param sample_data default NULL; character or data.frame; matching table used.
#' @return sample information table.
check_sample_table <- function(sample_data = NULL){
		# read according to the input class
		if(inherits(sample_data, "character")){
			if(grepl("csv", sample_data)){
				sample_data <- read.csv(sample_data, row.names = 1, stringsAsFactors = FALSE)
			}else{
				sample_data <- read.delim(sample_data, row.names = 1, stringsAsFactors = FALSE)
			}			
		}else{
			if(! inherits(sample_data, "data.frame")){
				stop("The input sample_data has unknown format! Must be character or data.frame format!")
			}
		}
	# output new abundance table
	sample_data
}

#' Get the website for a 'MetaCyc' pathway name
#'
#' @param pathway default NULL; character vector; one or more MetaCyc pathway names.
#' @return character vector.
#' @examples
#' metacyc_pathway_website("FOLSYN-PWY")
#' @export
metacyc_pathway_website <- function(pathway = NULL){
	if(is.null(pathway)){
		stop("Please input the pathway name!")
	}else{
		web_prefix <- "https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object="
		site <- paste0(web_prefix, pathway)
		site
	}
}
