#' Merge multiple microtable objects into one.
#'
#' @description
#' Merge multiple \code{microtable} objects into one unified \code{microtable} object (e.g. combining 16S bacteria data with ITS fungi data).
#' This is especially useful when performing cross-domain analysis.
#' 
#' @param ... two or more \code{microtable} objects to be merged.
#' @param prefix default NULL; a character vector of prefixes added to the feature names of each object 
#'   to avoid name conflicts. The length must be the same as the number of objects. 
#'   For example, \code{c("Bac_", "Fun_")} for combining bacterial and fungal data. 
#'   If NULL, no prefix is added (the user must ensure feature names do not conflict).
#' @param sample_intersect default TRUE; whether to take the intersection of samples across all objects. 
#'   If TRUE, only samples existing in all objects are retained. 
#'   If FALSE, the union of samples is taken, and missing samples are filled with 0.
#' @param auto_tidy default TRUE; whether to automatically invoke the \code{tidy_dataset} function on the merged object.
#' @return a merged \code{microtable} object.
#' @examples
#' \donttest{
#' # Example 1: Merge two microtable objects with prefixes
#' library(microeco)
#' data(otu_table_16S)
#' data(taxonomy_table_16S)
#' data(sample_info_16S)
#'
#' # Create simulated 16S and ITS microtable objects
#' m1 <- microtable$new(otu_table = otu_table_16S, 
#'   sample_table = sample_info_16S, 
#'   tax_table = taxonomy_table_16S)
#' m2 <- microtable$new(otu_table = otu_table_16S, 
#'   sample_table = sample_info_16S, 
#'   tax_table = taxonomy_table_16S)
#'
#' # Merge with prefixes
#' merged <- merge_microtable(m1, m2, prefix = c("Bac_", "Fun_"))
#'
#' # Example 2: Merge without prefixes (user ensures no name conflicts)
#' merged2 <- merge_microtable(m1, m2)
#'
#' # Example 3: Merge with union of samples instead of intersection
#' merged3 <- merge_microtable(m1, m2, sample_intersect = FALSE, 
#'   prefix = c("Bac_", "Fun_"))
#' }
#' @export
merge_microtable <- function(..., prefix = NULL, sample_intersect = TRUE, auto_tidy = TRUE){
	obj_list <- list(...)

	if(length(obj_list) < 2){
		stop("At least two microtable objects must be provided!")
	}

	obj_class <- unlist(lapply(obj_list, inherits, what = "microtable"))
	if(!all(obj_class)){
		stop("All input objects must be of class 'microtable'!")
	}

	has_otu <- unlist(lapply(obj_list, function(x) !is.null(x$otu_table) && nrow(x$otu_table) > 0))
	if(!all(has_otu)){
		stop("All input microtable objects must have a non-empty otu_table!")
	}

	if(!is.null(prefix)){
		if(length(prefix) != length(obj_list)){
			stop("The length of prefix must be the same as the number of objects! Provided: ", length(prefix), " but need: ", length(obj_list))
		}
	}

	n_objects <- length(obj_list)
	message("Merging ", n_objects, " microtable objects ...")

	for(i in seq_len(n_objects)){
		obj <- obj_list[[i]]

		if(!is.null(prefix)){
			old_otu_names <- rownames(obj$otu_table)
			new_otu_names <- paste0(prefix[i], old_otu_names)
			rownames(obj$otu_table) <- new_otu_names

			if(!is.null(obj$tax_table)){
				tax_keep <- rownames(obj$tax_table) %in% old_otu_names
				if(!all(tax_keep)){
					message("Removing ", sum(!tax_keep), " features from tax_table of object ", i, " that are not in otu_table ...")
					obj$tax_table <- obj$tax_table[tax_keep, , drop = FALSE]
				}
				otu_in_tax <- old_otu_names %in% rownames(obj$tax_table)
				if(!all(otu_in_tax)){
					stop("Some features in otu_table of object ", i, " are not found in tax_table!")
				}
				match_idx <- match(old_otu_names, rownames(obj$tax_table))
				rownames(obj$tax_table) <- new_otu_names[match_idx]
			}

			if(!is.null(obj$phylo_tree)){
				phylo_tips <- obj$phylo_tree$tip.label
				mt <- match(old_otu_names, phylo_tips)
				phylo_tips[mt[!is.na(mt)]] <- new_otu_names[!is.na(mt)]
				obj$phylo_tree$tip.label <- phylo_tips
			}

			if(!is.null(obj$rep_fasta)){
				fasta_names <- names(obj$rep_fasta)
				if(!is.null(fasta_names)){
					match_fa <- match(old_otu_names, fasta_names)
					if(any(is.na(match_fa))){
						warning("Some features in otu_table of object ", i, " are not found in rep_fasta names!")
					}
					fasta_names[match_fa[!is.na(match_fa)]] <- new_otu_names[!is.na(match_fa)]
					names(obj$rep_fasta) <- fasta_names
				}
			}

			message("Added prefix '", prefix[i], "' to features in object ", i, " ...")
		}

		obj_list[[i]] <- obj
	}

	if(sample_intersect){
		sample_sets <- lapply(obj_list, function(x) colnames(x$otu_table))
		common_samples <- Reduce(intersect, sample_sets)

		if(length(common_samples) == 0){
			stop("No common samples found across all objects! Please check the sample names. ",
				"If you want to keep all samples, set sample_intersect = FALSE.")
		}

		all_samples_unique <- unique(unlist(sample_sets))
		if(length(common_samples) < length(all_samples_unique)){
			removed_samples <- setdiff(all_samples_unique, common_samples)
			message("Using sample intersection mode. ", length(common_samples), 
				" common samples retained; ", length(removed_samples), 
				" samples removed: ", paste(removed_samples, collapse = ", "))
		}else{
			message("All ", length(common_samples), " samples are common across all objects.")
		}

		for(i in seq_len(n_objects)){
			obj_list[[i]]$otu_table <- obj_list[[i]]$otu_table[, common_samples, drop = FALSE]
		}
	}else{
		all_samples <- unique(unlist(lapply(obj_list, function(x) colnames(x$otu_table))))
		common_samples <- all_samples
		message("Using sample union mode. ", length(all_samples), " total samples retained.")

		for(i in seq_len(n_objects)){
			obj <- obj_list[[i]]
			existing_cols <- colnames(obj$otu_table)
			missing_cols <- setdiff(all_samples, existing_cols)
			if(length(missing_cols) > 0){
				fill_mat <- matrix(0, nrow = nrow(obj$otu_table), ncol = length(missing_cols))
				colnames(fill_mat) <- missing_cols
				obj$otu_table <- cbind(obj$otu_table, fill_mat)
				obj$otu_table <- obj$otu_table[, all_samples, drop = FALSE]
			}else{
				obj$otu_table <- obj$otu_table[, all_samples, drop = FALSE]
			}
			obj_list[[i]] <- obj
		}
	}

	merged_otu <- do.call(rbind, lapply(obj_list, function(x) x$otu_table))
	merged_otu <- as.data.frame(merged_otu)
	message("Merged otu_table: ", nrow(merged_otu), " features, ", ncol(merged_otu), " samples.")

	tax_list <- lapply(obj_list, function(x) x$tax_table)
	tax_list <- tax_list[!unlist(lapply(tax_list, is.null))]

	if(length(tax_list) > 0){
		all_tax_cols <- unique(unlist(lapply(tax_list, colnames)))
		for(i in seq_along(tax_list)){
			missing_cols <- setdiff(all_tax_cols, colnames(tax_list[[i]]))
			if(length(missing_cols) > 0){
				for(col in missing_cols){
					tax_list[[i]][[col]] <- "unclassified"
				}
			}
			tax_list[[i]] <- tax_list[[i]][, all_tax_cols, drop = FALSE]
		}
		merged_tax <- do.call(rbind, tax_list)
		merged_tax <- as.data.frame(merged_tax, stringsAsFactors = FALSE)
		message("Merged tax_table: ", nrow(merged_tax), " features, ", ncol(merged_tax), " taxonomic ranks.")
	}else{
		merged_tax <- NULL
		message("No tax_table in any input object; set to NULL.")
	}

	phylo_exists <- unlist(lapply(obj_list, function(x) !is.null(x$phylo_tree)))
	phylo_count <- sum(phylo_exists)

	if(phylo_count > 1){
		warning("Multiple phylo_tree objects found. Phylogenetic trees cannot be directly merged. Setting phylo_tree to NULL.")
		merged_phylo <- NULL
	}else if(phylo_count == 1){
		merged_phylo <- obj_list[[which(phylo_exists)[1]]]$phylo_tree
		message("Retaining phylo_tree from object ", which(phylo_exists)[1], ".")
	}else{
		merged_phylo <- NULL
	}

	fasta_list <- lapply(obj_list, function(x) x$rep_fasta)
	fasta_list <- fasta_list[!unlist(lapply(fasta_list, is.null))]

	if(length(fasta_list) > 1){
		fasta_classes <- unlist(lapply(fasta_list, class))
		fasta_types <- unlist(lapply(fasta_list, function(x) class(x)[1]))
		if(length(unique(fasta_types)) > 1){
			warning("rep_fasta objects have different types. Cannot merge. Setting rep_fasta to NULL.")
			merged_fasta <- NULL
		}else{
			merged_fasta <- do.call(c, fasta_list)
			message("Merged rep_fasta: ", length(merged_fasta), " sequences.")
		}
	}else if(length(fasta_list) == 1){
		merged_fasta <- fasta_list[[1]]
		message("Retaining rep_fasta from one object: ", length(merged_fasta), " sequences.")
	}else{
		merged_fasta <- NULL
	}

	base_sample <- obj_list[[1]]$sample_table
	if(!is.null(base_sample)){
		base_sample$..rownames.. <- rownames(base_sample)
		for(i in seq_len(n_objects)[-1]){
			if(!is.null(obj_list[[i]]$sample_table)){
				next_sample <- obj_list[[i]]$sample_table
				next_sample$..rownames.. <- rownames(next_sample)
				new_cols <- setdiff(colnames(next_sample), colnames(base_sample))
				if(length(new_cols) > 0){
					base_sample <- merge(base_sample, 
						next_sample[, c("..rownames..", new_cols), drop = FALSE], 
						by = "..rownames..", all.x = TRUE, sort = FALSE)
				}
			}
		}
		rownames(base_sample) <- base_sample$..rownames..
		base_sample$..rownames.. <- NULL
		base_sample <- base_sample[common_samples, , drop = FALSE]
		message("Merged sample_table: ", nrow(base_sample), " samples, ", ncol(base_sample), " columns.")
	}else{
		base_sample <- NULL
	}

	for(i in seq_len(n_objects)){
		message("Object ", i, ": ", nrow(obj_list[[i]]$otu_table), " features.")
	}

	newmt <- microtable$new(
		otu_table = merged_otu,
		sample_table = base_sample,
		tax_table = merged_tax,
		phylo_tree = merged_phylo,
		rep_fasta = merged_fasta
	)

	newmt$taxa_abund <- NULL
	newmt$alpha_diversity <- NULL
	newmt$beta_diversity <- NULL
	newmt$auto_tidy <- auto_tidy

	if(auto_tidy){
		newmt$tidy_dataset()
	}

	message("Merge completed! Total features: ", nrow(newmt$otu_table), 
		", total samples: ", nrow(newmt$sample_table))

	newmt
}
