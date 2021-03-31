# file2meco
Conversion between files from some tools and microtable object in microeco.


## Install file2meco

Install file2meco from github.

```r
# If devtools package is not installed, first install it
install.packages("devtools")
# then install microeco; https://github.com/ChiLiubio/microeco
devtools::install_github("ChiLiubio/microeco")
# then install file2meco
devtools::install_github("ChiLiubio/file2meco")
```

If the installation from github is failed because of the bad internet, download the packages first, then install them locally.

```r
devtools::install_local("microeco-master.zip")
devtools::install_local("file2meco-master.zip")
```

# Let's begin

```r
library(file2meco)
```


## QIIME files to microtable

The qiime1meco() function can be used to construct the microtable object using the raw OTU file from QIIME.

```r
# First install the dependence -- package qiimer
install.packages(system.file("extdata", "biom_0.3.12.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "qiimer_0.9.4.tar.gz", package="microeco"), repos = NULL, type = "source")
```

```r
# see the help document
?qiime1meco
# Let's run the examples
# use the raw data files stored inside the package
otu_file_path <- system.file("extdata", "otu_table_raw.txt", package="file2meco")
sample_file_path <- system.file("extdata", "sample_info.csv", package="file2meco")
phylo_file_path <- system.file("extdata", "rep_phylo.tre", package="file2meco")
# if you want to use Tax4Fun2 approach, you need read the representative sequences and add it to the microtable object.
rep_fasta_path <- system.file("extdata", "rep.fna", package="file2meco")
# contruct microtable object
qiime1meco(otu_table = otu_file_path, commented = FALSE, sample_data = sample_file_path)
qiime1meco(otu_table = otu_file_path, commented = FALSE, sample_data = sample_file_path, phylo_tree = phylo_file_path)
qiime1meco(otu_table = otu_file_path, commented = FALSE, sample_data = sample_file_path, phylo_tree = phylo_file_path, rep_fasta = rep_fasta_path)
```



## QIIME2 files to microtable

The qiime2meco() function can be used to construct the microtable object using files from QIIME2.
You need first install the required qiime2R package from github, see https://github.com/jbisanz/qiime2R


```r
# see the help document
?qiime2meco
# Let's run the examples
# use data files inside the package which were downloaded from (https://docs.qiime2.org/2020.8/tutorials/pd-mice/).
abund_file_path <- system.file("extdata", "dada2_table.qza", package="file2meco")
sample_file_path <- system.file("extdata", "sample-metadata.tsv", package="file2meco")
taxonomy_file_path <- system.file("extdata", "taxonomy.qza", package="file2meco")
phylo_file_path <- system.file("extdata", "tree.qza", package="file2meco")
rep_fasta_path <- system.file("extdata", "dada2_rep_set.qza", package="file2meco")
# contruct microtable object
qiime2meco(ASV_data = abund_file_path, sample_data = sample_file_path, taxonomy_data = taxonomy_file_path)
qiime2meco(ASV_data = abund_file_path, sample_data = sample_file_path, taxonomy_data = taxonomy_file_path, phylo_tree = phylo_file_path, rep_fasta = rep_fasta_path)
```


## Conversion between phyloseq and microtable
We provide two functions meco2phyloseq() and phyloseq2meco() for the conversion between microtable object and phyloseq object in phyloseq package.

```r
# Please first install phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
library(phyloseq)
```

```r
# from microtable to phyloseq object
data("dataset")
physeq <- meco2phyloseq(dataset)
physeq
```

```r
# from phyloseq to microtable object
data("GlobalPatterns")
meco_dataset <- phyloseq2meco(GlobalPatterns)
meco_dataset
```


## HUMAnN output files to microtable

The humann2meco() function can be used to construct the microtable object using metagenomic analysis files from HUMAnN2 and HUMAnN3.
Currently, it only support the KEGG pathway abundance file input. More input format will be supported.

```r
?humann2meco
# use the raw data files stored inside the package
abund_file_path <- system.file("extdata", "example_HUMAnN_KEGG_abund.tsv", package="file2meco")
sample_file_path <- system.file("extdata", "example_HUMAnN_sample_info.tsv", package="file2meco")
match_file_path <- system.file("extdata", "example_HUMAnN_match_table.tsv", package="file2meco")
humann2meco(abund_table = abund_file_path)
humann2meco(abund_table = abund_file_path, sample_data = sample_file_path, match_table = match_file_path)
```

```r
# Let's try more interesting usages with microeco
library(file2meco)
library(microeco)
library(magrittr)
test <- humann2meco(abund_table = abund_file_path, sample_data = sample_file_path, match_table = match_file_path)
test$tax_table %<>% subset(level1 != "unclassified")
test$tidy_dataset()
# rel = FALSE donot use relative abundance
test$cal_abund(select_cols = 1:3, rel = FALSE)
test1 <- trans_abund$new(test, taxrank = "level2", ntaxa = 10)
test1$plot_bar(facet = "Group", ylab_title = "Abundance (RPK)")
# select both function and taxa
test$cal_abund(select_cols = c("level1", "Phylum", "Genus"), rel = TRUE)
test1 <- trans_abund$new(test, taxrank = "Phylum", ntaxa = 10, delete_part_prefix = T)
test1$plot_bar(facet = "Group")
# functional biomarker
test$cal_abund(select_cols = 1:3, rel = TRUE)
test1 <- trans_diff$new(test, method = "lefse", group = "Group")
test1$plot_lefse_bar(LDA_score = 3)
# taxa biomarker
test$cal_abund(select_cols = 4:9, rel = TRUE)
test1 <- trans_diff$new(test, method = "lefse", group = "Group")
test1$plot_lefse_bar(LDA_score = 2)
```

# Other tools

This package will be updated continuously......  
Any idea/suggestion will be considered. We also appreciate that anyone can join us to make the package better.


## References
  - Chi Liu, Yaoming Cui, Xiangzhen Li, Minjie Yao, microeco: an R package for data mining in microbial community ecology, FEMS Microbiology Ecology, Volume 97, Issue 2, February 2021, fiaa255.
  - McMurdie PJ, Holmes S (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS ONE 8(4): e61217. 
  - Bolyen, E., Rideout, J.R., Dillon, M.R. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852–857 (2019).
