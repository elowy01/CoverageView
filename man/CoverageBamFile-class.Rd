\name{CoverageBamFile-class}
\Rdversion{1.1}
\docType{class}
\alias{CoverageBamFile-class}
\alias{CoverageBamFile}
% accessors
\alias{path,CoverageBamFile-method}
\alias{reads_mapped,CoverageBamFile-method}
\alias{run_type,CoverageBamFile-method}

\alias{show,CoverageBamFile-method}

\title{Class \code{"CoverageBamFile"}}
\description{
The \code{CoverageBamFile} class contains information on a BAM file and inherits fields from the \code{\linkS4class{BamFile}} class in the \code{Rsamtools} package.
}

\section{Objects from the Class}{
Use \code{CoverageBamFile()} to create a reference to a BAM file.

This class represents an alignment in the BAM format that can be processed by the different methods of the CoverageView package.

The \code{reads_mapped} argument is optional and will be automatically calculated (if it is not provided) when the \code{normalization} argument is set in the \code{cov.matrix} and \code{cov.interval} functions
}
\arguments{
       
    \item{\code{path}}{ A character string that details the path to the BAM file }
    \item{\code{reads_mapped}}{ A number representing the number of reads having an alignment. Default = 0 }
    \item{\code{run_type}}{ A character string describing the type of sequencing run. Possible values are 'single' and 'paired'. Default='single' }
    \item{...}{Additional arguments}
}

\references{
BAM format specification: http://samtools.sourceforge.net/SAMv1.pdf
}
\author{
Ernesto Lowy <ernestolowy@gmail.com>
}

\seealso{
  \itemize{
      \item \code{\link{BamFile}} 
      \item \code{\link{CoverageBigWigFile}}
      }
}
\examples{
#get a Bam test file
inputfile<-system.file("extdata","treat.bam",package="CoverageView")
  
#create a CoverageBamFile object
trm<-CoverageBamFile(inputfile)
  
#create a CoverageBamFile object including the information of the number
#of reads aligned and the sequencing run type
trm1<-CoverageBamFile(inputfile,reads_mapped=28654321,run_type="single")
}
\keyword{classes}
