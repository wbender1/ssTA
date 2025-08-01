# Shape Space Trait Association Algorithm

![R]
![GitHub](https://img.shields.io/badge/license-MIT-green)

An algorithm for determining association between viral-gene coding-sequences and categorical traits, such as mild/severe disease. 
Utilized for peer-reviewed publication: https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-024-01930-7
Also includes functions for performing Antigenic Cartography (Sequences-Based Mapping).  
https://pmc.ncbi.nlm.nih.gov/articles/PMC5809904/

## ✨ Features
- **Quality Control**: Cleans sequences to desired Nucleic Acid lengths
- **Translation**: Performs translation of Nuceic Acids into Amino Acids
- **Alignment**: Aligns Amino Acids
- **Antigenic Hamming Distances**: Calculates the Antigenic Hamming Distances between sequences and produces a distance matrix
- **Principal Coordinate Analysis (PCoA)**: Performs a dimension reduction on distance matrix to produce lower-dimensional data
- **Trait Association**: Determines association between categorical trait and Antigenic Hamming Distance matrix.  
Produces a:
- Stats file containing anosim, adonis2, bDis.perm, and bDis.anova pvalues.
- Pdf file containing an Anosim plot.
- Pdf file containing a bDisper sdEllipse.





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
