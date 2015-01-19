\name{cov.matrix}
\alias{cov.matrix}
\alias{cov.matrix,ANY,missing-method}
\alias{cov.matrix,ANY,ANY-method}
\title{Compute a coverage matrix}
\description{
This method generates a coverage matrix for which each column is a different genomic feature (i.e. the coordinate of the transcription start site) or genomic interval (defined by start/end coordinates) and each row can be either a position or a bin related to this feature or interval. Each element in the matrix is the result of computing the coverage value for this position/bin that is at a certain distance upstream or dowstream the genomic feature (i.e. the TSS) or it is in the genomic interval
}
\usage{
\S4method{cov.matrix}{ANY,missing}(tr,coordfile,normalization,bin_width,extend,no_windows,offset,num_cores)
\S4method{cov.matrix}{ANY,ANY}(tr,ctl,coordfile,normalization,bin_width,extend,no_windows,offset,do,num_cores)
}
\arguments{
\item{tr}{An instance of a \code{\linkS4class{CoverageBamFile}} or \code{\linkS4class{CoverageBigWigFile}} object}
\item{ctl}{An instance of a \code{\linkS4class{CoverageBamFile}} or \code{\linkS4class{CoverageBigWigFile}} object. This argument is optional and if provided then this object will be used in the arithmetic operations with the coverages specified with the \code{do} argument. Default=NULL}
\item{coordfile}{File in the \code{BED} or \code{GFF} format containing the coordinates used to generate the coverage matrix}
\item{normalization}{Normalize the coverage values in the matrix. Possible normalization types are \code{rpm} (Read per million). Default=NULL}
\item{bin_width}{Calculate the average coverage for bins of length defined by this argument. Default is 1 nucleotide}
\item{extend}{Controls the number of nucleotides to extend relative to the anchor coordinates (i.e. TSS) as specified in the \code{coordfile} argument}
\item{no_windows}{Number of windows used to divide the interval defined by the start and end coordinates as specified in the \code{coordfile} argument}
\item{offset}{If the \code{no_windows} argument is set, this controls the number of windows to extend the interval defined by the start/end coordinates. Default=0}
\item{do}{Used to set the arithmetic operation to perform in a certain position or bin. Possible values are: \code{log2ratio}, \code{ratio}, \code{substraction}. This argument only makes sense when either a \code{\linkS4class{CoverageBamFile}} or \code{\linkS4class{CoverageBigWigFile}} object is provided with the \code{ctl} argument. Default= \code{log2ratio}}
\item{num_cores}{Number of cores to use. Default= The \code{detectCores()} function from the library \code{parallel} will be used to check the number of cores available in the system and this number will be used as the default value}
}
\details{
The \code{cov.matrix} function receives one or two \code{\linkS4class{CoverageBamFile}} or \code{\linkS4class{CoverageBigWigFile}} objects and returns a coverage matrix in which every column can be either a single genomic feature (i.e. the coordinates of the transcription start site for a gene) or a certain genomic interval defined by the start and end coordinates specified in the file used with the \code{coordfile} argument, and each row represents a different position/bin relative to the feature or a certain window within a genomic interval. This behaviour is controlled by the \code{extend} and \code{no_windows} parameters. For example, if one set the \code{extend} argument, then each row represents a single genomic position (or bin with a a size defined by the \code{bin_width} argument) relative to the anchor coordinate. Each matrix element will be the result of calculating the depth of coverage for that particular position/bin relative to the anchor coordinate defined in \code{coordfile}. If the \code{bin_width} parameter is >1 nt, then the average coverate for each bin is used for all downstream calculations. On the other hand, if the \code{no_windows} parameter is set, then each row in the matrix will represent a certain genomic window within an interval defined by the start and end coordinates as specified in the \code{coordfile} file. In this case, each element in the matrix will be the average coverage value for that particular genomic window. The number of windows in which the genomic interval is divided is controlled by the \code{no_windows} argument. Finally, if only one \code{\linkS4class{CoverageBamFile}} or \code{\linkS4class{CoverageBigWigFile}} is provided using the \code{tr} argument, then the matrix elements are the result of computing the coverage value or the average coverage value for a certain position or a certain genomic bin/window respectively. If, on the other hand, two files are provided using the \code{tr} and the \code{ctl} arguments, then the matrix elements will be the result of performing a certain arithmetic operation (specified with the \code{do} argument) with the coverages.

If the \code{normalization} parameter is set to 'rpm' then the number of reads aligned into the reference is required to perform the normalization. This number can be provided using the \code{reads_mapped} argument of the \code{CoverageBamFile} or \code{CoverageBigWigFile} objects, or it is automatically calculated if the function is used with a \code{CoverageBamFile} object. NOTE: If the \code{reads_mapped} argument is set then the automatic calculation will NOT be performed and the provided number will be taken as correct for all downstream calculations
}
\value{
Returns a \code{matrix}
}
\references{
BAM format specification: http://samtools.sourceforge.net/SAMv1.pdf

WIG format specification: https://genome.ucsc.edu/FAQ/FAQformat.html

BED format specification: https://genome.ucsc.edu/FAQ/FAQformat.html
}
\author{
Ernesto Lowy <ernestolowy@gmail.com>
}
\seealso{
\code{\link{draw.profile}}
\code{\link{draw.heatmap}}
}
\examples{
###########
##BAM files
###########
  
## Generating a coverage matrix with a single BAM file
  
#get a BAM test file
aBAMfile<-system.file("extdata","treat.bam",package="CoverageView")
  
#get a BED file with the TSS (transcription start site) coordinates of some genes
bedfile<-system.file("extdata","testTSS.bed",package="CoverageView")
  
#create the CoverageBamFile object
trm<-CoverageBamFile(aBAMfile,reads_mapped=28564510)
  
#generate the coverage matrix extending 100 nts on each side of the provided TSS
#in the bedfile
DF1=cov.matrix(tr=trm,coordfile=bedfile,extend=100,num_cores=2)
  
#generate the coverage matrix extending 100 nts on each side of the TSS provided
#in the bedfile and normalize the resulting coverages
DF2=cov.matrix(tr=trm,coordfile=bedfile,extend=100,normalization="rpm",num_cores=2)
  
#generate the coverage matrix extending 100 nts on each side of the TSS provided
#in the bedfile and normalize the resulting coverages. This time we calculate the average
#coverage in bins of 2 nucleotides
DF3=cov.matrix(tr=trm,coordfile=bedfile,extend=100,normalization="rpm",bin_width=2,num_cores=2)
  
## Generating a coverage matrix with 2 BAM files
  
#get 2 BAM test files
treatBAMfile<-system.file("extdata","treat.bam",package="CoverageView")
ctrlBAMfile<-system.file("extdata","ctrl.bam",package="CoverageView")
  
#get a BED file with the TSS (transcription start site) coordinates of some genes
bedfile<-system.file("extdata","testTSS.bed",package="CoverageView")
  
#create 2 CoverageBamFile objects
trm<-CoverageBamFile(treatBAMfile,reads_mapped=28564510)
ctl<-CoverageBamFile(ctrlBAMfile,reads_mapped=26713667)
  
#generate the coverage matrix extending 100 nts on each side of the TSS provided
#in the bedfile and normalize the resulting coverages.The matrix elements are obtained
#by computing the ratio of the coverages of the trm against the ctl objects and then
#calculating the log2 of the ratios
DF4=cov.matrix(tr=trm,ctl=ctl,coordfile=bedfile,extend=100,normalization="rpm",do="log2ratio",num_cores=2)
  
#####################
  
#get a treatment BAM test file
treatBAMfile<-system.file("extdata","treat.bam",package="CoverageView")
  
#get a GFF file with the chr,start and end coordinates of different genomic
#features (i.e. CDS)
gffile<-system.file("extdata","test.gff",package="CoverageView")
  
#create the 'treatment' CoverageBamFile object
trm<-CoverageBamFile(treatBAMfile,reads_mapped=28564510)
  
#generate the coverage matrix dividing each genomic interval defined by the start and
#end coordinates in the gff file into 10 windows and calculating the average coverage
#for each window
DF1=cov.matrix(trm,coordfile=gffile,no_windows=10,num_cores=2)
  
#generate the coverage matrix dividing each genomic interval defined by the start and
#end coordinates in the gff file into 10 windows and calculating the average coverage
#for each window, this time we extend the genomic interval by 1 window before the start
#coordinate and 1 window after the end coordinate (offset argument is set to 1)
DF1=cov.matrix(trm,coordfile=gffile,no_windows=10,offset=1)
  
###########
##BigWIG files
###########
  
## Generating a coverage matrix with a single WIG file
  
#get a bigWIG test file
abigWIGfile<-system.file("extdata","treat.bw",package="CoverageView")

#get a BED file with the TSS (transcription start site) coordinates of some genes
bedfile<-system.file("extdata","testTSS.bed",package="CoverageView")

#create the CoverageBigWigFile object
trm<-CoverageBigWigFile(abigWIGfile)
  
#generate the coverage matrix extending 100 nts on each side of the TSS provided
DF1=cov.matrix(trm,coordfile=bedfile,extend=100,bin_width=10,num_cores=2)
  
## Generating a coverage matrix with 2 BigWIG files
  
#get 2 BigWIG test files
treatBigWIGfile<-system.file("extdata","treat.bw",package="CoverageView")
ctrlBigWIGfile<-system.file("extdata","ctrl.bw",package="CoverageView")
  
#create 2 CoverageBigWigFile objects
trm<-CoverageBigWigFile(treatBigWIGfile)
ctl<-CoverageBigWigFile(ctrlBigWIGfile)
  
#generate the coverage matrix extending 100 nts on each side of the TSS provided
#in the bedfile .The matrix elements are obtained by computing the ratio of the 
#coverages of the trm against the ctl objects and then calculating the log2 of the ratios
DF2=cov.matrix(tr=trm,ctl=ctl,coordfile=bedfile,extend=100,do="log2ratio",num_cores=2)
}

