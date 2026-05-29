
# file2meco <a href="https://chiliubio.github.io/microeco_tutorial/"><img src="https://user-images.githubusercontent.com/20815519/128602544-78d53642-b445-4686-a22a-1ef3c0726ce7.png" width=150 align="right" ></a>

Convert files of some tools to microtable object of microeco package


[![CRAN](https://www.r-pkg.org/badges/version/file2meco)](https://cran.r-project.org/web/packages/file2meco/index.html)
[![CRAN](https://cranlogs.r-pkg.org/badges/grand-total/file2meco)](https://cran.r-project.org/web/packages/file2meco/index.html)
![](https://img.shields.io/badge/Release-v0.9.1-blue.svg) ![](https://img.shields.io/badge/Test-v1.0.0-red.svg)


## Install file2meco

Install file2meco package from CRAN.

```r
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages("file2meco", repos = BiocManager::repositories())
```

## Citation

Chi Liu, Xiangzhen Li, Felipe R. P. Mansoldo, Tong Chen, Fanzheng Meng, Ruixiang Tang, Siyu Zhou, Qinghua Yang, Ruixin Shao & Minjie Yao. 
microeco 2: A comprehensive R package for downstream analysis of microbiome omics data. 
**iMeta**, 2026, 5: e70132. https://doi.org/10.1002/imt2.70132

Chi Liu, Felipe R. P. Mansoldo, Hankang Li, Alane Beatriz Vermelho, Raymond Jianxiong Zeng, Xiangzhen Li & Minjie Yao. 
A workflow for statistical analysis and visualization of microbiome omics data using the R microeco package. 
**Nature Protocols**, 2026, 21: 1300–1324. https://doi.org/10.1038/s41596-025-01239-4


## Tutorial

https://chiliubio.github.io/microeco_tutorial/file2meco-package.html


## Supported tools

Currently supported tools for data conversion include:  
QIIME, QIIME2, HUMAnN, MetaPhlAn, Kraken2/Bracken, Ncyc, PICRUSt2, Tax4Fun/Tax4Fun2, ViromeScan, R package phyloseq and TreeSummarizedExperiment.

## Contributing

We welcome any contribution \! 
Any idea/suggestion will be considered.
By participating in this project you agree to abide by the terms outlined in the [Contributor Code of Conduct](CONDUCT.md).



## References
  - https://chiliubio.github.io/microeco_tutorial/references.html#references
