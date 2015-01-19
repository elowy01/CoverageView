\name{CoverageBigWigFile-class}
\Rdversion{1.1}
\docType{class}
\alias{CoverageBigWigFile-class}
\alias{CoverageBigWigFile}

\alias{show,CoverageBigWigFile-method}

\title{Class \code{"CoverageBigWigFile"}}

\description{
The CoverageBigWigFile class contains information on a BigWIG file and inherits fields from the \code{\linkS4class{BigWigFile}} class in the \code{rtracklayer} package.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{CoverageBigWigFile()}. This class represents an alignment in the BigWIG format for which each genomic position has a depth of coverage value associated.

The \code{reads_mapped} argument is optional but it is required if the \code{normalization} parameter is set in the \code{cov.matrix} and \code{cov.interval} functions

}
\section{Slots}{
  \describe{ 
    \item{\code{resource}}{A character string that details the path to the BigWIG file }
    \item{\code{reads_mapped}}{A number representing the number of reads aligned in the file. Default = 0}
    \item{...}{Additional arguments}
  }
}
\references{
WIG format specification: https://genome.ucsc.edu/FAQ/FAQformat.html
}


\author{
Ernesto Lowy <ernestolowy@gmail.com>
}

\seealso{
 \itemize{
      \item \code{\link{BigWigFile}} 
      \item \code{\link{CoverageBamFile}}
      }
}
\examples{
##get BigWIG test file
inputfile<-system.file("extdata","treat.bw",package="CoverageView")

#create a CoverageBigWigFile object
trm<-CoverageBigWigFile(inputfile)

#create a CoverageBigWigFile object including the information of the number of reads aligned
trm1<-CoverageBigWigFile(inputfile,reads_mapped=28654321)

##get BigWIG test file
inputfile<-system.file("extdata","treat.bw",package="CoverageView")

#create a CoverageBigWigFile object
trm<-CoverageBigWigFile(inputfile)

#create a CoverageBigWigFile object including the information of the number of reads aligned
trm1<-CoverageBigWigFile(inputfile,reads_mapped=28654321)
}

