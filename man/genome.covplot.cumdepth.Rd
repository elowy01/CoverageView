\name{genome.covplot.cumdepth}
\alias{genome.covplot.cumdepth}
\alias{genome.covplot.cumdepth,CoverageBamFile-method}
\alias{genome.covplot.cumdepth,list-method}
\title{Generates a cumulative genome coverage plot}
\description{
  This method generates a plot showing the percentage of the genome covered
at different read depths
}
\usage{
  \S4method{genome.covplot.cumdepth}{CoverageBamFile}(data,outfile,max_depth)
  \S4method{genome.covplot.cumdepth}{list}(data,outfile,max_depth)
}
\arguments{
   \item{data}{Either an instance of \code{\linkS4class{CoverageBamFile}} or a list of \code{\linkS4class{CoverageBamFile}} objects}
  \item{outfile}{URL of the \code{.png} file where the plot will be created}
  \item{max_depth}{Maximum read depth to be displayed in the X-axis}
}
\details{
This method receives either a single \code{\linkS4class{CoverageBamFile}} object or a list of \code{CoverageBamFile} objects and generates a plot for which the X-axis is a range of cumulative read depths and the Y-axis is the percentage of the genome covered at a certain read depth. If a list of \code{CoverageBamFile} objects is passed to the function then it will generate a different coloured line for each of the passed objects
}
\author{
  Ernesto Lowy <ernestolowy@gmail.com>
}
\seealso{
  \code{\link{genome.covplot.depth}}
}
\examples{
##draw a cumulative coverage plot for a test case BAM file
  
#get a BAM test file
treatBAMfile<-system.file("extdata","treat.bam",package="CoverageView")
  
#create a CoverageBamFile object
trm<-CoverageBamFile(treatBAMfile)
  
#draw the plot
genome.covplot.cumdepth(trm,outfile="test.png")
  
#draw the plot setting the max_depth parameter (30X in this case)
genome.covplot.cumdepth(trm,outfile="test.png",max_depth=30)
  
##draw two overlapping cumulative coverage plots for two different BAM files
  
#get the first BAM file
treatBAMfile<-system.file("extdata","treat.bam",package="CoverageView")
#create the CoverageBamFile object
trm<-CoverageBamFile(treatBAMfile)
  
#get the second BAM test file
ctrlBAMfile<-system.file("extdata","ctrl.bam",package="CoverageView")
#create the CoverageBamFile object
ctl<-CoverageBamFile(ctrlBAMfile)

#create a list with the two files
input_d=list(trm,ctl)
  
#draw the plot
genome.covplot.cumdepth(input_d,outfile="test.png")
}
  