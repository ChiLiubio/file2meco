#' Transform 'QIIME2' results to 'microtable' object.
#'
#' @description
#' Transform 'QIIME2' qza results to microtable object.
#' @param ASV_data the ASV data, such as the \code{'data2_table.qza'}.
#' @param sample_data default NULL; the sample metadata table; four types of formats are available: 
#'   1) q2-type tab seperated file of QIIME2, such as the \code{'sample-metadata.tsv'} in the example;
#'   2) comma seperated file with the suffix csv or tab seperated file with suffix tsv/txt;
#'   3) Excel type file with the suffix xlsx or xls; require \code{readxl} package to be installed;
#'   4) \code{data.frame} object from R.
#' @param taxonomy_data default NULL; the taxonomy data, such as the \code{'taxonomy.qza'}.
#' @param phylo_tree default NULL; the phylogenetic tree, such as the \code{'tree.qza'}.
#' @param rep_fasta default NULL; the representative sequences, such as the \code{'dada2_rep_set.qza'}.
#' @param ... parameter passed to \code{microtable$new} function of \code{microeco} package, such as \code{auto_tidy} parameter.
#' @return \code{microtable} object.
#' @examples
#' \dontrun{
#' # The data files is downloaded from https://docs.qiime2.org/2020.8/tutorials/pd-mice/ 
#' #   and stored inside the package.
#' abund_file_path <- system.file("extdata", "dada2_table.qza", package="file2meco")
#' sample_file_path <- system.file("extdata", "sample-metadata.tsv", package="file2meco")
#' taxonomy_file_path <- system.file("extdata", "taxonomy.qza", package="file2meco")
#' qiime2meco(ASV_data = abund_file_path)
#' qiime2meco(ASV_data = abund_file_path, sample_data = sample_file_path)
#' qiime2meco(ASV_data = abund_file_path, sample_data = sample_file_path, 
#'   taxonomy_data = taxonomy_file_path)
#' }
#' @export
qiime2meco <- function(ASV_data, sample_data = NULL, taxonomy_data = NULL, phylo_tree = NULL, rep_fasta = NULL, ...){
	# if(!require(qiime2R)){
		# stop("qiime2R package not installed!")
	# }
	# Read ASV data
	ASV <- as.data.frame(read_qza(ASV_data)$data)
	#  Read metadata
	if(!is.null(sample_data)){
		if(is_q2metadata(sample_data)){
			sample_data <- read_q2metadata(sample_data)
			rownames(sample_data) <- as.character(sample_data[, 1])
		}else{
			sample_data <- check_sample_table(sample_data = sample_data)
		}
	}
	# Read taxonomy table
	if(!is.null(taxonomy_data)){
		taxonomy_data <- read_qza(taxonomy_data)
		taxonomy_data <- q2_parse_taxonomy(taxonomy_data$data)
		# Make the taxonomic table clean, this is very important.
		taxonomy_data %<>% tidy_taxonomy
	}

	# Read phylo tree
	if(!is.null(phylo_tree)){
		phylo_tree <- read_qza(phylo_tree)$data
	}
	if(!is.null(rep_fasta)){
		rep_fasta <- read_qza(rep_fasta)$data
	}
	dataset <- microtable$new(sample_table = sample_data, tax_table = taxonomy_data, otu_table = ASV, phylo_tree = phylo_tree, rep_fasta = rep_fasta, ...)
	dataset
}



####################################################################
# The following functions come from the 'qiime2R' package https://github.com/jbisanz/qiime2R
# All rights reserved.
# The reason implementing those from qiime2R is that submitting package to CRAN have a strict requirement on the repository.

###############################################################

# read 'qiime2' metadata (.tsv)
# 
# Loads a 'qiime2' metadata file wherein the 2nd line contains the #q2:types line dictating the type of variable (categorical/numeric)
#
# param file path to the input file, ex: file="~/data/moving_pictures/table.qza"

# return a data.frame wherein the first column is SampleID
read_q2metadata <- function(file) {
  defline<-suppressWarnings(readLines(file)[2])
  defline<-strsplit(defline, split="\t")[[1]]
  
  defline[grep("numeric", tolower(defline))]<-"double"
  defline[grep("categorical|q2:types", tolower(defline))]<-"factor"
  defline[defline==""]<-"factor"
  
  coltitles<-strsplit(suppressWarnings(readLines(file)[1]), split='\t')[[1]]
  metadata<-read.table(file, header=F, col.names=coltitles, skip=2, sep='\t', colClasses = defline, check.names = FALSE)
  colnames(metadata)[1]<-"SampleID"
  
  return(metadata)
}

# checks if metadata is in 'qiime2' (.tsv)
#
# Checks to see if a file is in 'qiime2' metadata format, ie contains #q2:types line dictating the type of variable (categorical/numeric)
#
# param file path to the input file, ex: file="~/data/moving_pictures/table.qza"

# return TRUE/FALSE
is_q2metadata <- function(file){
  suppressWarnings(
  if(grepl("^#q2:types", readLines(file)[2])){return(TRUE)}else{return(FALSE)}
  )
}


# read 'qiime2' artifacts (.qza)
#
# extracts embedded data and object metadata into an R session
#
# file path to the input file, ex: file="~/data/moving_pictures/table.qza"
# tmp a temporary directory that the object will be decompressed to (default="tempdir()")
# rm should the decompressed object be removed at completion of function (T/F default=TRUE)
# return a named list of objects.
read_qza <- function(file, tmp, rm) {

if(missing(tmp)){tmp <- tempdir()}
if(missing(file)){stop("Path to artifact (.qza) not provided.")}
if(!file.exists(file)){stop("Input artifact (",file,") not found. Please check path and/or use list.files() to see files in current working directory.")}
if(missing(rm)){rm=TRUE} #remove the decompressed object from tmp
if(!grepl("qza$", file)){stop("Provided file is not qiime2 artifact (.qza).")}

unzip(file, exdir=tmp)
unpacked<-unzip(file, exdir=tmp, list=TRUE)

artifact<-read_yaml(paste0(tmp,"/", paste0(gsub("/..+","", unpacked$Name[1]),"/metadata.yaml"))) #start by loading in the metadata not assuming it will be first file listed
artifact$contents<-data.frame(files=unpacked)
artifact$contents$size=sapply(paste0(tmp, "/", artifact$contents$files), file.size)
artifact$version=read.table(paste0(tmp,"/",artifact$uuid, "/VERSION"))

  #get data dependent on format
  #get data dependent on format
if(grepl("BIOMV", artifact$format)){
  artifact$data<-read_q2biom(paste0(tmp, "/", artifact$uui,"/data/feature-table.biom"))
} else if (artifact$format=="NewickDirectoryFormat"){
  artifact$data<-read.tree(paste0(tmp,"/",artifact$uuid,"/data/tree.nwk"))
} else if (artifact$format=="DistanceMatrixDirectoryFormat") {
  artifact$data<-as.dist(read.table(paste0(tmp,"/", artifact$uuid, "/data/distance-matrix.tsv"), header=TRUE, row.names=1, fill= TRUE))
} else if (grepl("StatsDirFmt", artifact$format)) {
  if(paste0(artifact$uuid, "/data/stats.csv") %in% artifact$contents$files.Name){artifact$data<-read.csv(paste0(tmp,"/", artifact$uuid, "/data/stats.csv"), header=TRUE, row.names=1)}
  if(paste0(artifact$uuid, "/data/stats.tsv") %in% artifact$contents$files.Name){artifact$data<-read.table(paste0(tmp,"/", artifact$uuid, "/data/stats.tsv"), header=TRUE, row.names=1, sep='\t')} #can be tsv or csv
} else if (artifact$format=="TSVTaxonomyDirectoryFormat"){
  artifact$data<- read.table(paste0(tmp,"/", artifact$uuid, "/data/taxonomy.tsv"), sep='\t', header=TRUE, quote="", comment.char="")
} else if (artifact$format=="DNASequencesDirectoryFormat") {
  artifact$data<- Biostrings::readDNAStringSet(paste0(tmp,"/",artifact$uuid,"/data/dna-sequences.fasta"))
} else if (artifact$format=="AlignedDNASequencesDirectoryFormat") {
  artifact$data<- Biostrings::readDNAMultipleAlignment(paste0(tmp,"/",artifact$uuid,"/data/aligned-dna-sequences.fasta"))
} else if (grepl("EMPPairedEndDirFmt|EMPSingleEndDirFmt|FastqGzFormat|MultiplexedPairedEndBarcodeInSequenceDirFmt|MultiplexedSingleEndBarcodeInSequenceDirFmt|PairedDNASequencesDirectoryFormat|SingleLanePerSamplePairedEndFastqDirFmt|SingleLanePerSampleSingleEndFastqDirFmt", artifact$format)) {
  artifact$data<-data.frame(files=list.files(paste0(tmp,"/", artifact$uuid,"/data")))
  artifact$data$size<-format(sapply(artifact$data$files, function(x){file.size(paste0(tmp,"/",artifact$uuid,"/data/",x))}, simplify = TRUE))
} else {
  message("Format not supported, only a list of internal files and provenance is being imported.")
  artifact$data<-list.files(paste0(tmp,"/",artifact$uuid, "/data"))
}

#Add Provenance
pfiles<-paste0(tmp,"/", grep("..+provenance/..+action.yaml", unpacked$Name, value=TRUE))
artifact$provenance<-lapply(pfiles, read_yaml)
names(artifact$provenance)<-grep("..+provenance/..+action.yaml", unpacked$Name, value=TRUE)

if(rm==TRUE){unlink(paste0(tmp,"/", artifact$uuid), recursive=TRUE)}

return(artifact)
}


# read 'qiime2' biom file (version 2.1)
#
# Loads a version 2.1 spec biom file (http://biom-format.org/documentation/format_versions/biom-2.1.html) as expected to be found within a 'qiime2' artifact.
#
# file path to the input file, ex: file="~/Downloads/3372d9e0-3f1c-43d8-838b-35c7ad6dac89/data/feature-table.biom"

# return a matrix of values
read_q2biom <- function(file) {
  if(missing(file)){stop("Path to biom file given")}
  if(!file.exists(file)){stop("File not found")}
  
  hdata<-h5read(file,"/")
  
  ftable<-
    sparseMatrix(
      p=hdata$observation$matrix$indptr,
      j=hdata$observation$matrix$indices,
      x=as.numeric(hdata$observation$matrix$data),
      index1=FALSE,
      dims=c(length(hdata$observation$ids), length(hdata$sample$ids)),
      dimnames=list(hdata$observation$ids,hdata$sample$ids)
    )
  
  return(as.matrix(ftable))
}

# Parse 'Qiime2' taxonomy
#
# taxonomy a table-like object containing the columns Feature.ID and Taxon. Can be imported using read_qza(file)$data.
# tax_sep The separator between taxonomic levels. Defaults to one compatible with both GreenGenes and SILVA ("; " OR ";")
# trim_extra Remove leading characters from taxonomic levels: ex: k__ or D_0__. TRUE/FALSE. default=TRUE 
# 
# Note: Assumes an assignment has been made to all levels. Fills missing assignments with NA.
# return a data.frame with feature IDs as row names and the columns: Kingdom, Phylum, Class, Order, Family, Genus, Species
#
q2_parse_taxonomy <- function(taxonomy, tax_sep, trim_extra){
  if(missing(taxonomy)){stop("Taxonomy Table not supplied.")}
  if(missing(trim_extra)){trim_extra=TRUE}
  if(missing(tax_sep)){tax_sep="; |;"}
  if(sum(colnames(taxonomy) %in% c("Feature.ID","Taxon"))!=2){stop("Table does not match expected format. ie does not have columns Feature.ID and Taxon.")}

  taxonomy<-taxonomy[,c("Feature.ID","Taxon")]
  if(trim_extra){
  taxonomy$Taxon<-gsub("[kpcofgs]__","", taxonomy$Taxon) #remove leading characters from GG
  taxonomy$Taxon<-gsub("D_\\d__","", taxonomy$Taxon) #remove leading characters from SILVA
  }
  taxonomy<-suppressWarnings(taxonomy %>% separate(Taxon, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep=tax_sep))
  taxonomy<-apply(taxonomy, 2, function(x) if_else(x=="", NA_character_, x))
  taxonomy<-as.data.frame(taxonomy)
  rownames(taxonomy)<-taxonomy$Feature.ID
  taxonomy$Feature.ID<-NULL
  return(taxonomy)  
}
