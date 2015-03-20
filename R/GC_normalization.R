library(Rsamtools)
library(CoverageView)
library(BSgenome)

# 1) read in .BAM file
bf <- BamFile("/Users/ernesto/test/SRR1737292_1.aln.sorted.bam")

reads <- GenomicAlignments::readGAlignmentsFromBam(bf)

lengths<-seqlengths(reads)

cov<-coverage(reads)

# work only with chr1
chr1_len<-lengths[1]

# 2) BIN the genome: create bins of 100000 in this case
intervals<-seq(from=0,to=chr1_len,100000)

#split intervals into each 2-element chunks
max<-2
# x will be a vector containing a seq that goes from 1 to length of intervals
x <- seq_along(intervals)
#list composed of each interval
d1 <- split(intervals, ceiling(x/max))

#get last element of list
last<-tail(d1,n=1)

if(length(last[[1]]==1)) {
    vlast<-unlist(last)
    vlast<-append(vlast,chr1_len)
    d1[[length(d1)]]<-vlast
}

# 3) Calculate mean coverage for each bin
rleO <- cov[["chr1"]]

calcMeanCov <- function(l,rleO) {
    #l is the list with intervals and rleO is object with coverages
    stretch=rleO[l[1]:l[2]]
    mean(stretch)
}
# returns a list with mean coverages per bin
res <- lapply(d1, calcMeanCov,rleO)

# 4) Calculate GC content per window
hg19 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
chr1<-hg19$chr1

calcGC <- function(l,chr) {
    browser()
    stretch=chr[l[1]:l[2]]
    counts<-alphabetFrequency(stretch)
}

res1 <- lapply(d1, calcGC,chr1)

