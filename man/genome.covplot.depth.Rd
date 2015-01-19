\name{genome.covplot.depth}
\alias{genome.covplot.depth}
\alias{genome.covplot.depth,CoverageBamFile-method}
\alias{genome.covplot.depth,list-method}
\title{
Generates a genome coverage plot
}
\description{
This method generates a plot showing the number of genomic positions reaching a certain read depth
}
\usage{
  \S4method{genome.covplot.depth}{CoverageBamFile}(data,outfile,max_depth)
  \S4method{genome.covplot.depth}{list}(data,outfile,max_depth)
}
\arguments{
  \item{data}{Either an instance of a \code{\linkS4class{CoverageBamFile}} object or a list of \code{\linkS4class{CoverageBamFile}} objects}
  \item{outfile}{URL of the \code{.png} file where the plot will be created}
  \item{max_depth}{Maximum read depth to be displayed in the X-axis}
}
\details{
This method receives either a single \code{\linkS4class{CoverageBamFile}} object or a list of \code{CoverageBamFile} objects and generates a plot for which the X-axis represents a range of coverage read depths and the Y-axis corresponds to the number of megabases having a specific read coverage value. If a list of \code{CoverageBamFile} objects is passed to the function then it will generate a different coloured line for each of the objects
}
\author{
  Ernesto Lowy <ernestolowy@gmail.com>
}
\seealso{
  \code{\link{genome.covplot.cumdepth}}
}
\examples{
##draw a coverage plot for a test case BAM file
  
#get a BAM test file
treatBAMfile<-system.file("extdata","treat.bam",package="CoverageView")
  
#create the CoverageBamFile object
trm<-CoverageBamFile(treatBAMfile)
  
#draw the plot
genome.covplot.depth(trm,outfile="test.png")
  
#draw the plot setting the max_depth parameter (70X in this case)
genome.covplot.depth(trm,outfile="test.png",max_depth=70)
  
##draw two overlapping coverage plots for two different test BAM files
  
#get a first BAM test file
treatBAMfile<-system.file("extdata","treat.bam",package="CoverageView")
#create the CoverageBamFile object
trm<-CoverageBamFile(treatBAMfile)
  
#get a second BAM test file
ctrlBAMfile<-system.file("extdata","ctrl.bam",package="CoverageView")
#create the CoverageBamFile object
ctl<-CoverageBamFile(ctrlBAMfile)

#create a list with the two files
input_d=list(trm,ctl)
  
#draw the plot
genome.covplot.depth(input_d,outfile="test.png")
}
