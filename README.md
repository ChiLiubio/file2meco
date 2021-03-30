# file2meco
Conversion between files from some tools and microtable.


## Install file2meco

Install file2meco from github.

```r
# If devtools package is not installed, first install it
install.packages("devtools")
# then install microeco
devtools::install_github("ChiLiubio/microeco")
# then install file2meco
devtools::install_github("ChiLiubio/file2meco")
```

If the installation from github is failed because of the bad internet, download the packages first, then install it locally.

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
# The data files is downloaded from https://docs.qiime2.org/2020.8/tutorials/pd-mice/ and stored inside the package.
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

## References
  - Chi Liu, Yaoming Cui, Xiangzhen Li, Minjie Yao, microeco: an R package for data mining in microbial community ecology, FEMS Microbiology Ecology, Volume 97, Issue 2, February 2021, fiaa255.
  - McMurdie PJ, Holmes S (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS ONE 8(4): e61217. 



