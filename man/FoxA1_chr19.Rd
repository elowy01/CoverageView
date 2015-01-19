\name{FoxA1_chr19}
\alias{FoxA1_chr19}
\docType{data}
\title{
Example of a coverage matrix using the ChIP-seq data for the FoxA1 transcription factor experiment
}
\description{
A coverage matrix calculated using the FoxA1 ChIP-seq data from the T-47D cell line. This matrix was computed using only the reads aligned into the Human chromosome 19
}
\usage{data(FoxA1_chr19)}
\format{
	\code{matrix} with the coverages for the genomic interval corresponding to the gene body of a certain set of genes
}
\details{
This matrix was obtained by running the \code{CoverageView} function named \code{cov.matrix} and setting the \code{extend} parameter to 1000 and the \code{bin\_width} paremeter to 10. The BED format coordinate file was passed to the function through the \code{coordfile} argument and contained the transcription start sites (TSS) of the genes that had a valid MACS peak in the region spanning 2500 nucleotides on each side of the TSS. All peaks, regardless its false discovery rate, were considered. For illustrative purposes only the genes located in the chr19 were analyzed.
}

\references{
Zhang et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol (2008) vol. 9 (9) pp. R137".
}

\source{
FoxA1 ChIP-seq data set for the T-47D cell line from the HudsonAlpha Institute that was used as a test case in the Nature protocol named:
Identifying ChIP-seq enrichment using MACS. Feng J, Liu T, Qin B, Zhang Y, Liu XS. Nat Protoc. 2012 Sep;7(9):1728-40.
}
\keyword{datasets}
