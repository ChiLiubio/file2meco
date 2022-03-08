
# file2meco <a href="https://chiliubio.github.io/microeco_tutorial/"><img src="https://user-images.githubusercontent.com/20815519/128602544-78d53642-b445-4686-a22a-1ef3c0726ce7.png" width=150 align="right" ></a>

Convert files of some tools to microtable object of microeco package


[![CRAN](https://www.r-pkg.org/badges/version/file2meco)](https://cran.r-project.org/web/packages/file2meco/index.html)
[![CRAN](https://cranlogs.r-pkg.org/badges/grand-total/file2meco)](https://cran.r-project.org/web/packages/file2meco/index.html)
![](https://img.shields.io/badge/Release-v0.2.0-blue.svg) ![](https://img.shields.io/badge/Test-v0.2.2-red.svg)


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

## Tutorial

Please see https://chiliubio.github.io/microeco_tutorial/file2meco-package.html for the examples or the help documents.

### QIIME files to microtable

The qiime1meco() function can be used to construct the microtable object using the raw OTU file from QIIME.

### QIIME2 files to microtable

The qiime2meco() function can be used to create the microtable object using files from QIIME2.

### HUMAnN metagenomic results to microtable

HUMAnN is an excellent tool for functional profiling analysis of metagenomes and metatranscriptomes at species-level (https://doi.org/10.1038/s41592-018-0176-y).
The humann2meco() function can be used to creat the microtable object using metagenomic analysis files from HUMAnN2 and HUMAnN3 (https://huttenhower.sph.harvard.edu/humann).
Currently, it supports both the MetaCyc and KEGG pathway abundance file input.

### MetaPhlAn
MetaPhlAn is an software used for metagenomic taxonomic profiling (https://doi.org/10.1038/nmeth.3589).
The format of MetaPhlAn classification results is usually called 'mpa' format.
The mpa2meco function is developed for this format conversion to microtable object.
See the following example of Kraken2 part.


### Kraken2
Kraken is a taxonomic sequence classifier that assigns taxonomic labels to DNA sequences.
Kraken examines the k-mers within a query sequence and uses the information within those k-mers to query a database. 
That database maps k-mers to the lowest common ancestor (LCA) of all genomes known to contain a given k-mer.
Kraken2 is the newest version (https://doi.org/10.1186/s13059-019-1891-0).
The merged Kraken2 results can be obtained by merge_metaphlan_tables.py from MetaPhlAn or combine_mpa.py from KrakenTools (https://ccb.jhu.edu/software/krakentools/).

### Ncyc

Ncyc database is a curated integrative database for fast and accurate metagenomic profiling of nitrogen cycling genes (https://doi.org/10.1093/bioinformatics/bty741).
The ncyc2meco() function is designed for construct the microtable object using gene abundance files from Ncyc (https://github.com/qichao1984/NCyc).


### Conversion between phyloseq and microtable
Two functions meco2phyloseq() and phyloseq2meco() were provided for the conversion between microtable object and phyloseq object in phyloseq package.

## Other tools

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
