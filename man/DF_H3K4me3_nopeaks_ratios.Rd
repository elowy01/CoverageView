\name{DF_H3K4me3_nopeaks_ratios}
\alias{DF_H3K4me3_nopeaks_ratios}
\docType{data}
\title{
Example of a matrix with the ratio of the coverages using the ChIP-seq data for the H3K4me3 histone modification experiment
}
\description{
A matrix with the ratio of the coverages calculated using the H3K4me3 ChIP-seq data from the K562 cell line and its respective control replicate. This matrix was computed using only the reads aligned into the Human chromosome 19
}
\usage{data(DF_H3K4me3_nopeaks_ratios)}
\format{
	\code{matrix} with the normalized (read per million) coverages for the genomic interval corresponding to the gene body of a certain set of genes
}
\details{
This matrix was obtained by running the \code{CoverageView} function named \code{cov.matrix} with two CoverageBamFile objects (one sample and one control respectively) and setting the \code{no_windows} parameter to 100, the \code{normalization} parameter to "rpm" and the offset parameter to 10. The \code{do} argument was set to \code{ratio}, so each element in the matrix will be the result of calculating the ratio of the coverages between the sample and the control.
The BED format coordinate file was passed to the function through the \code{coordfile} argument and contained the genomic coordinates of a set of genes for which the region starting 2500 nucleotides before the transcription start site (TSS) and ending in the transcription end site (TES) did not contain any enrichment peak identified by MACS using the ChIP-seq data. For illustrative purposes only the genes located in the chr19 were analyzed.
}

\references{
Zhang et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol (2008) vol. 9 (9) pp. R137".
}

\source{
H3K4me3 ChIP-seq data set for the K562 cell line from the HudsonAlpha Institute that was used as a test case in the Nature protocol named:
Identifying ChIP-seq enrichment using MACS. Feng J, Liu T, Qin B, Zhang Y, Liu XS. Nat Protoc. 2012 Sep;7(9):1728-40.
}
\keyword{datasets}
