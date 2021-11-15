
# file2meco <a href="https://chiliubio.github.io/microeco_tutorial/"><img src="https://user-images.githubusercontent.com/20815519/128602544-78d53642-b445-4686-a22a-1ef3c0726ce7.png" width=150 align="right" ></a>

Convert files of some tools to microtable object of microeco package


[![CRAN](https://www.r-pkg.org/badges/version/file2meco)](https://cran.r-project.org/web/packages/file2meco/index.html)
[![CRAN](https://cranlogs.r-pkg.org/badges/grand-total/file2meco)](https://cran.r-project.org/web/packages/file2meco/index.html)
![](https://img.shields.io/badge/Release-v0.2.0-blue.svg) ![](https://img.shields.io/badge/Test-v0.2.1-red.svg)


## Install file2meco

Install file2meco package from CRAN directly.

```r
install.packages("file2meco")
```

Or install the latest development version from github.

```r
# If devtools package is not installed, first install it
install.packages("devtools")
devtools::install_github("ChiLiubio/file2meco")
```

# Let's begin

```r
library(file2meco)
```


## QIIME files to microtable

The qiime1meco() function can be used to construct the microtable object using the raw OTU file from QIIME.


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
qiime1meco(otu_table = otu_file_path, commented = FALSE)
qiime1meco(otu_table = otu_file_path, commented = FALSE, sample_data = sample_file_path)
qiime1meco(otu_table = otu_file_path, commented = FALSE, sample_data = sample_file_path, phylo_tree = phylo_file_path)
qiime1meco(otu_table = otu_file_path, commented = FALSE, sample_data = sample_file_path, phylo_tree = phylo_file_path, rep_fasta = rep_fasta_path)
```



## QIIME2 files to microtable

The qiime2meco() function can be used to construct the microtable object using files from QIIME2.


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
# construct microtable object
qiime2meco(ASV_data = abund_file_path)
qiime2meco(ASV_data = abund_file_path, sample_data = sample_file_path, taxonomy_data = taxonomy_file_path)
qiime2meco(ASV_data = abund_file_path, sample_data = sample_file_path, taxonomy_data = taxonomy_file_path, phylo_tree = phylo_file_path, rep_fasta = rep_fasta_path)
```


## HUMAnN metagenomic results to microtable

HUMAnN is an excellent tool for functional profiling analysis of metagenomes and metatranscriptomes at species-level (https://doi.org/10.1038/s41592-018-0176-y).
The humann2meco() function can be used to creat the microtable object using metagenomic analysis files from HUMAnN2 and HUMAnN3 (https://huttenhower.sph.harvard.edu/humann).
Currently, it supports both the MetaCyc and KEGG pathway abundance file input.


```r
library(file2meco)
library(microeco)
library(magrittr)
?humann2meco
sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", package="file2meco")
match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")

# MetaCyc pathway database based analysis
# use the raw data files stored inside the package for MetaCyc pathway database based analysis
abund_file_path <- system.file("extdata", "example_HUMAnN_MetaCyc_abund.tsv", package="file2meco")
# the default db is "MetaCyc"
humann2meco(abund_table = abund_file_path, db = "MetaCyc")
humann2meco(abund_table = abund_file_path, db = "MetaCyc", sample_data = sample_file_path, match_table = match_file_path)
# Let's try more interesting usages with microeco
test <- humann2meco(abund_table = abund_file_path, db = "MetaCyc", sample_data = sample_file_path, match_table = match_file_path)
test$tidy_dataset()
# rel = FALSE donot use relative abundance
test$cal_abund(select_cols = 1:3, rel = FALSE)
test$taxa_abund$Superclass1 %<>% .[!grepl("unclass", rownames(.)), ]
test1 <- trans_abund$new(test, taxrank = "Superclass1", ntaxa = 10)
test1$plot_bar(facet = "Group", ylab_title = "Abundance (RPK)")
# select both function and taxa
test$cal_abund(select_cols = c("Superclass1", "Phylum", "Genus"), rel = TRUE)
test1 <- trans_abund$new(test, taxrank = "Phylum", ntaxa = 10, delete_part_prefix = T)
test1$plot_bar(facet = "Group")
# functional biomarker
test$cal_abund(select_cols = 1:3, rel = TRUE)
test$taxa_abund$Superclass1 %<>% .[!grepl("unclass", rownames(.)), ]
test1 <- trans_diff$new(test, method = "lefse", group = "Group")
test1$plot_lefse_bar(use_number = 1:20)
# taxa biomarker
test$cal_abund(select_cols = 4:9, rel = TRUE)
test$taxa_abund$Phylum %<>% .[!grepl("unclass", rownames(.)), ]
test1 <- trans_diff$new(test, method = "lefse", group = "Group")
test1$plot_lefse_bar(LDA_score = 2)
```


```r
# use KEGG pathway based HUMAnN result
abund_file_path <- system.file("extdata", "example_HUMAnN_KEGG_abund.tsv", package="file2meco")
test <- humann2meco(abund_table = abund_file_path, db = "KEGG", sample_data = sample_file_path, match_table = match_file_path)
test$tax_table %<>% subset(level1 != "unclassified")
test$tidy_dataset()
# rel = FALSE donot use relative abundance
test$cal_abund(select_cols = 1:3, rel = FALSE)
test1 <- trans_abund$new(test, taxrank = "level2", ntaxa = 10)
test1$plot_bar(facet = "Group", ylab_title = "Abundance (RPK)")
# select both function and taxa
test$cal_abund(select_cols = c("level1", "Phylum", "Genus"), rel = TRUE)
test$taxa_abund$level1 %<>% .[!grepl("unclass", rownames(.)), ]
test$taxa_abund$Phylum %<>% .[!grepl("unclass", rownames(.)), ]
test1 <- trans_abund$new(test, taxrank = "Phylum", ntaxa = 10, delete_part_prefix = T)
test1$plot_bar(facet = "Group")
# functional biomarker
test$cal_abund(select_cols = 1:3, rel = TRUE)
test1 <- trans_diff$new(test, method = "lefse", group = "Group")
test1$plot_lefse_bar(LDA_score = 3)
```


## MetaPhlAn
MetaPhlAn is an software used for metagenomic taxonomic profiling (https://doi.org/10.1038/nmeth.3589).
The format of MetaPhlAn classification results is usually called 'mpa' format.
The mpa2meco function is developed for this format conversion to microtable object.
See the following example of Kraken2 part.


## Kraken2
Kraken is a taxonomic sequence classifier that assigns taxonomic labels to DNA sequences.
Kraken examines the k-mers within a query sequence and uses the information within those k-mers to query a database. 
That database maps k-mers to the lowest common ancestor (LCA) of all genomes known to contain a given k-mer.
Kraken2 is the newest version (https://doi.org/10.1186/s13059-019-1891-0).
The merged Kraken2 results can be obtained by merge_metaphlan_tables.py from MetaPhlAn or combine_mpa.py from KrakenTools (https://ccb.jhu.edu/software/krakentools/).

```r
# the example is metagenomic classification result
# use the raw data files stored inside the package
abund_file_path <- system.file("extdata", "example_kraken2_merge.txt", package="file2meco")
sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", 
  package="file2meco")
match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
library(microeco)
library(file2meco)
library(magrittr)
mpa2meco(abund_table = abund_file_path)
test <- mpa2meco(abund_table = abund_file_path, sample_data = sample_file_path, 
  match_table = match_file_path)
test$tidy_dataset()
```


## Ncyc

Ncyc database is a curated integrative database for fast and accurate metagenomic profiling of nitrogen cycling genes (https://doi.org/10.1093/bioinformatics/bty741).
The ncyc2meco() function is designed for construct the microtable object using gene abundance files from Ncyc (https://github.com/qichao1984/NCyc).


```r
?ncyc2meco
# use the raw data files stored inside the package
abund_file_path <- system.file("extdata", "example_Ncyc_table.tsv", package="file2meco")
sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", package="file2meco")
match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
ncyc2meco(abund_table = abund_file_path)
ncyc2meco(abund_table = abund_file_path, sample_data = sample_file_path, match_table = match_file_path)
```

```r
# Let's try more interesting usages with microeco
library(file2meco)
library(microeco)
library(magrittr)
test <- ncyc2meco(abund_table = abund_file_path, sample_data = sample_file_path, match_table = match_file_path)
test$tidy_dataset()
# use split_group = TRUE to calculate the pathway abundance with multipe map correspondance
test$cal_abund(select_cols = 1:2, rel = TRUE, split_group = TRUE, split_column = "Pathway")
test$taxa_abund$Pathway %<>% .[!grepl("unclass", rownames(.)), ]
test1 <- trans_abund$new(test, taxrank = "Pathway")
test1$plot_bar(bar_type = "notfull")
# for gene abundance, no splitting on the pathways
test$cal_abund(select_cols = 1:2, rel = TRUE, split_group = FALSE)
test$taxa_abund$Gene %<>% .[!grepl("unclass", rownames(.)), ]
test1 <- trans_abund$new(test, taxrank = "Gene")
test1$plot_bar(bar_type = "notfull")
```

## Conversion between phyloseq and microtable
Two functions meco2phyloseq() and phyloseq2meco() were provided for the conversion between microtable object and phyloseq object in phyloseq package.

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


# Other tools

Other converting methods will be developed continuously......  

## Contributing

We welcome any contribution \! 
Any idea/suggestion will be considered.
By participating in this project you agree to abide by the terms outlined in the [Contributor Code of Conduct](CONDUCT.md).



## References
  - Chi Liu, Yaoming Cui, Xiangzhen Li, Minjie Yao, microeco: an R package for data mining in microbial community ecology, FEMS Microbiology Ecology, Volume 97, Issue 2, February 2021, fiaa255.
  - Bolyen, E., Rideout, J.R., Dillon, M.R. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852–857 (2019).
  - Franzosa EA, McIver LJ, Rahnavard G, Thompson LR, Schirmer M, Weingart G, Schwarzberg Lipson K, Knight R, Caporaso JG, Segata N, Huttenhower C. Species-level functional profiling of metagenomes and metatranscriptomes. Nat Methods 15: 962-968 (2018).
  - Truong, D., Franzosa, E., Tickle, T. et al. MetaPhlAn2 for enhanced metagenomic taxonomic profiling. Nat Methods 12, 902–903 (2015). https://doi.org/10.1038/nmeth.3589
  - Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0
  - Qichao Tu, Lu Lin, Lei Cheng, Ye Deng, Zhili He, NCycDB: a curated integrative database for fast and accurate metagenomic profiling of nitrogen cycling genes, Bioinformatics, Volume 35, Issue 6, 15 March 2019, 1040–1048
  - Caspi et al 2020, "The MetaCyc database of metabolic pathways and enzymes - a 2019 update", Nucleic Acids Research 48(D1):D445-D453
  - Segata, N., Izard, J., Waldron, L., Gevers, D., Miropolsky, L., Garrett, W. S., & Huttenhower, C. (2011). Metagenomic biomarker discovery and explanation. Genome Biology, 12(6), R60.
  - McMurdie PJ, Holmes S (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS ONE 8(4): e61217. 
  - https://www.genome.jp/kegg/
