#' Get the website for a MetaCyc pathway name
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
