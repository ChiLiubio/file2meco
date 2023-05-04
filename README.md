
# file2meco <a href="https://chiliubio.github.io/microeco_tutorial/"><img src="https://user-images.githubusercontent.com/20815519/128602544-78d53642-b445-4686-a22a-1ef3c0726ce7.png" width=150 align="right" ></a>

Convert files of some tools to microtable object of microeco package


[![CRAN](https://www.r-pkg.org/badges/version/file2meco)](https://cran.r-project.org/web/packages/file2meco/index.html)
[![CRAN](https://cranlogs.r-pkg.org/badges/grand-total/file2meco)](https://cran.r-project.org/web/packages/file2meco/index.html)
![](https://img.shields.io/badge/Release-v0.6.0-blue.svg) ![](https://img.shields.io/badge/Test-v0.6.1-red.svg)


## Install file2meco

Install file2meco package from CRAN.

```r
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages("file2meco", repos = BiocManager::repositories())
```

## Citation
Liu, C., Li, X., Mansoldo, F.R.P., An, J., Kou, Y., Zhang, X., Wang, J., Zeng, J., Vermelho, A.B., Yao, M., 2022. 
Microbial habitat specificity largely affects microbial co-occurrence patterns and functional profiles in wetland soils. 
Geoderma 418, 115866. https://doi.org/10.1016/j.geoderma.2022.115866


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
  - https://chiliubio.github.io/microeco_tutorial/references.html#references
