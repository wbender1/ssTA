# ssTA

To install the package locally:
1. install.packages("devtools")
2. library(devtools)
3. devtools::install_github("wbender1/ssTA")

To install required packages if necessary:
- install.packages("seqinr")
- install.packages("bioseq")
- install.packages("stringdist")
- install.packages("Rtsne")
- install.packages("vegan")
- if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
  BiocManager::install("msa")

To obtain tutorial directory:
git clone https://github.com/wbender1/ssTA.git
