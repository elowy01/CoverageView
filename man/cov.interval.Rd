\name{cov.interval}
\alias{cov.interval}
\alias{cov.interval,CoverageBamFile-method}
\alias{cov.interval,CoverageBigWigFile-method}
\alias{cov.interval,CoverageBamFile,CoverageBamFile-method}
\alias{cov.interval,CoverageBigWigFile,CoverageBigWigFile-method}
\title{arithmetic operation on an interval}
\description{
This method computes a numeric vector with the result of a certain arithmetic operation with the coverage for two particular BAM or BigWIG files within a specific genomic interval defined by the user
}
\usage{
  \S4method{cov.interval}{CoverageBamFile,CoverageBamFile}(tr,ctl,normalization,chr,start,end,bin_width,do)
  \S4method{cov.interval}{CoverageBigWigFile,CoverageBigWigFile}(tr,ctl,normalization,chr,start,end,bin_width,do)
}
\arguments{
  \item{tr}{An instance of a \code{\linkS4class{CoverageBamFile}} or a \code{\linkS4class{CoverageBigWigFile}} object used to compute an arithmetic operation with the coverage values for the genomic interval defined by the \code{chr}, \code{start} and \code{end} parameters}
  \item{ctl}{An instance of a \code{\linkS4class{CoverageBamFile}} or a \code{\linkS4class{CoverageBigWigFile}} object used to compute a certain arithmetic operation with the coverage values for the genomic interval defined by the \code{chr}, \code{start} and \code{end} parameters}
  \item{normalization}{Normalize the coverage values in the resulting numeric vector. Possible normalization types are \code{rpm} (Read per million). Default=NULL}
  \item{chr}{chromosome of the genomic interval}
  \item{start}{start coordinate of the genomic interval. If the \code{start} and \code{end} arguments are not set, then the entire chromosome defined with the \code{chr} argument will be used for the calculations}
  \item{end}{end coordinate of the genomic interval}
  \item{bin_width}{calculate the average coverage within bins of size defined by this argument. Default is 1 nucleotide}
  \item{do}{specify the arithmetic operation to perform on the genomic interval. Possible values are 'log2ratio','ratio' or 'substraction'.Default=\code{log2ratio}}
}
\details{
The \code{cov.interval} function receives 2 \code{CoverageBamFiles} or \code{CoverageBigWigFiles} objects and returns a numeric vector with the result of doing a certain arithmetic operation with the two files using the coverage values for a certain position or genomic bin that is included in the genomic interval defined by the user using the \code{chr}, \code{start} and \code{end} arguments. The size of the bins is controlled by the \code{bin_width} argument and for each bin the average coverage value is computed and used for the different arithmetic operations.The arithmetic operation to perform is set with the \code{do} argument.

If the \code{normalization} parameter is set to 'rpm' then the number of reads aligned into the reference is required to perform the normalization. This number can be provided using the \code{reads_mapped} argument of the \code{CoverageBamFile} or \code{CoverageBigWigFile} objects, or it is automatically calculated if the function is used with a \code{CoverageBamFile} object. NOTE: If the \code{reads_mapped} argument is set then the automatic calculation will NOT be performed and the provided number will be taken as correct for all downstream calculations.
}
\references{
BAM format specification: http://samtools.sourceforge.net/SAMv1.pdf

WIG format specification: https://genome.ucsc.edu/FAQ/FAQformat.html
}
\author{
Ernesto Lowy <ernestolowy@gmail.com>
}
\examples{
##BAM files
  
#get treatment and control test files
treat1file<-system.file("extdata","treat.bam",package="CoverageView")
control1file<-system.file("extdata","ctrl.bam",package="CoverageView")
  
#create two CoverageBamFile objects representing single-end alignments
trm1<-CoverageBamFile(treat1file,run_type="single")
ctl1<-CoverageBamFile(control1file,run_type="single")
  
#calculate the ratio of the coverages for the defined genomic interval using a bin_width equal to 10 nts
cov1=cov.interval(trm1,ctl1,chr="chrI",start=1,end=100,bin_width=10,do="ratio")
  
#create a WIG file with the obtained vector with the ratios
outfolder=system.file("extdata",package="CoverageView")
an_outfile1=paste(outfolder,"out.wig",sep="/")
export.wig(cov1,outfile=an_outfile1)
  
##BigWIG files
  
#get a treatment and control test files
treat2file<-system.file("extdata","treat.bw",package="CoverageView")
control2file<-system.file("extdata","ctrl.bw",package="CoverageView")
  
#create the 'treatment' and 'control' CoverageBigWigFile objects
trm2<-CoverageBigWigFile(path=treat2file,reads_mapped=28564510)
ctl2<-CoverageBigWigFile(path=control2file,reads_mapped=26713667)
  
#calculate the ratio of the coverages for the defined genomic interval
cov2=cov.interval(trm2,ctl2,bin_width=1,chr="chrI",start=1,end=1000,do='ratio')
}