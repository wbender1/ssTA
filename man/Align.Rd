% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Align.R
\name{Align}
\alias{Align}
\title{Alignment}
\usage{
Align(input_file, output_file)
}
\arguments{
\item{input_file}{name of a fasta file containing Amino Acid sequences}

\item{output_file}{name of the output fasta file}
}
\value{
a fasta file containing aligned sequences
}
\description{
Aligns amino acid sequences using msa ClustalOmega.

Required package: msa
}
\examples{
setwd("~/ssTA/tutorial/")
myAlignedSeqs <- Align("Output_GOOD/PracticeSeqs_Translated_GOOD.fasta", "Output_GOOD/PracticeSeqs_Aligned_GOOD.fasta")

}
