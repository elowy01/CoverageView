### R code from vignette source 'CoverageView.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=40)


###################################################
### code chunk number 2: preliminaries
###################################################
library(CoverageView)


###################################################
### code chunk number 3: example-CoverageBamFile
###################################################
treatBAMfile<-system.file("extdata","treat.bam",package="CoverageView")
trm<-CoverageBamFile(treatBAMfile,reads_mapped=28564510)


###################################################
### code chunk number 4: example-Wigfile
###################################################
treatBigWIGfile<-system.file("extdata","treat.bw",package="CoverageView")
trm<-CoverageBigWigFile(treatBigWIGfile,reads_mapped=28564510)


###################################################
### code chunk number 5: example-cov.matrix-1 (eval = FALSE)
###################################################
## # get the FoxA1 BAM file for chr19
## FoxA1_chr19_BAM_url="http://www.ebi.ac.uk/~ernesto/FoxA1.chr19.bam"
## download.file(FoxA1_chr19_BAM_url,"./FoxA1.chr19.bam")
## # get the FoxA1 BAM index file
## FoxA1_chr19_BAI_url="http://www.ebi.ac.uk/~ernesto/FoxA1.chr19.bam.bai"
## download.file(FoxA1_chr19_BAI_url,"./FoxA1.chr19.bam.bai")
## ##
## # now instantiate a CoverageBamFile object
## trm<-CoverageBamFile("./FoxA1.chr19.bam",reads_mapped=168042)
## tss_chr19_bed<-system.file("extdata","FoxA1.chr19.bed",package="CoverageView")


###################################################
### code chunk number 6: example-cov.matrix-1b
###################################################
data(FoxA1_chr19)


###################################################
### code chunk number 7: example-cov.matrix-2 (eval = FALSE)
###################################################
## FoxA1_chr19<-cov.matrix(trm,coordfile=tss_chr19_bed,extend=1000,num_cores=2,
## bin_width=10)


###################################################
### code chunk number 8: example-cov.matrix-2b
###################################################
FoxA1_chr19[1:3,1:5]


###################################################
### code chunk number 9: example-cov.matrix-profile
###################################################
draw.profile(FoxA1_chr19,ylab="avg coverage",outfile="FoxA1_chr19.png",
main="FoxA1 chr19")


###################################################
### code chunk number 10: example-cov.matrix-profile
###################################################
write.profile(FoxA1_chr19,outfile="FoxA1_chr19.txt")


###################################################
### code chunk number 11: example-cov.matrix-heatmap
###################################################
draw.heatmap(FoxA1_chr19,outfile="FoxA1_chr19_heatmap.png")


###################################################
### code chunk number 12: example-cov.matrix-histone-1 (eval = FALSE)
###################################################
## # get the sample BAM file for chr19
## H3K4me3_BAM_url="http://www.ebi.ac.uk/~ernesto/H3K4me3.chr19.bam"
## download.file(H3K4me3_BAM_url,"./H3K4me3.chr19.bam")
## # get also the index for the previous file
## H3K4me3_BAI_url="http://www.ebi.ac.uk/~ernesto/H3K4me3.chr19.bam.bai"
## download.file(H3K4me3_BAI_url,"./H3K4me3.chr19.bam.bai")
## 
## #get the control BAM file for chr19
## H3K4me3_Control_BAM_url="http://www.ebi.ac.uk/~ernesto/H3K4me3_Control.chr19.bam"
## download.file(H3K4me3_Control_BAM_url,"./H3K4me3_Control.chr19.bam")
## # get also the index for the previous file
## H3K4me3_Control_BAI_url="http://www.ebi.ac.uk/~ernesto/H3K4me3_Control.chr19.bam.bai"
## download.file(H3K4me3_Control_BAI_url,"./H3K4me3_Control.chr19.bam.bai")


###################################################
### code chunk number 13: example-cov.matrix-histone-2 (eval = FALSE)
###################################################
## trm<-CoverageBamFile("./H3K4me3.chr19.bam",reads_mapped=864924)
## ctl<-CoverageBamFile("./H3K4me3_Control.chr19.bam",reads_mapped=319369)


###################################################
### code chunk number 14: example-cov.matrix-histone-3
###################################################
H3K4me3_chr19_bed<-system.file("extdata","H3K4me3.chr19.bed",
package="CoverageView")


###################################################
### code chunk number 15: example-cov.matrix-histone-4 (eval = FALSE)
###################################################
## DF_H3K4me3<-cov.matrix(trm,coordfile=H3K4me3_chr19_bed,no_windows=100,
## offset=10,num_cores=2,normalization="rpm")
## DF_H3K4me3_ctl<-cov.matrix(ctl,coordfile=H3K4me3_chr19_bed,no_windows=100,
## offset=10,num_cores=2,normalization="rpm")


###################################################
### code chunk number 16: example-cov.matrix-histone-4b
###################################################
data(DF_H3K4me3)
data(DF_H3K4me3_ctl)


###################################################
### code chunk number 17: example-cov.matrix-histone-5
###################################################
input_list=list(DF_H3K4me3,DF_H3K4me3_ctl) 
draw.profile(data=input_list,ylab="avg coverage",outfile="H3K4me3cmp.png"
,main="H3K4me3")


###################################################
### code chunk number 18: example-cov.matrix-histone-6
###################################################
input_list=list(DF_H3K4me3,DF_H3K4me3_ctl) 
draw.heatmap(data=input_list,outfile="H3K4me3cmp_heatmap.png")


###################################################
### code chunk number 19: example-cov.matrix-histone-7
###################################################
data(DF_H3K36me3)
data(DF_H3K36me3_control)


###################################################
### code chunk number 20: example-cov.matrix-histone-8
###################################################
input_list=list(DF_H3K4me3,DF_H3K4me3_ctl,DF_H3K36me3,DF_H3K36me3_control) 
draw.heatmap(data=input_list,outfile="H3K4me3VSH3K36me3_heatmap.png")


###################################################
### code chunk number 21: example-cov.matrix-histone-9 (eval = FALSE)
###################################################
## H3K4me3_chr19_nopeaks_bed<-system.file("extdata","H3K4me3.chr19.nopeaks.bed",
## package="CoverageView")
## DF_H3K4me3_nopeaks_ratios=cov.matrix(trm,ctl,coordfile=H3K4me3_chr19_nopeaks_bed,
## no_windows=100,offset=10,num_cores=2,normalization="rpm",do="log2ratio")
## save(DF_H3K4me3_nopeaks_ratios,file="DF_H3K4me3_nopeaks_ratios.RData")


###################################################
### code chunk number 22: example-cov.matrix-histone-9b
###################################################
data(DF_H3K4me3_nopeaks_ratios)


###################################################
### code chunk number 23: example-cov.matrix-histone-10
###################################################
draw.profile(DF_H3K4me3_nopeaks_ratios,ylab="avg coverage",
outfile="H3K4me3nopeaks_ratios.png",main="H3K4me3 ratio")


###################################################
### code chunk number 24: example-cov.interval-histone-1 (eval = FALSE)
###################################################
## # get the sample BAM file for chr19
## H3K4me3_BAM_url="http://www.ebi.ac.uk/~ernesto/H3K4me3.chr19.bam"
## download.file(H3K4me3_BAM_url,"./H3K4me3.chr19.bam")
## # get also the index file for previous file
## H3K4me3_BAI_url="http://www.ebi.ac.uk/~ernesto/H3K4me3.chr19.bam.bai"
## download.file(H3K4me3_BAI_url,"./H3K4me3.chr19.bam.bai")
## 
## #get the control BAM file for chr19
## H3K4me3_Control_BAM_url="http://www.ebi.ac.uk/~ernesto/H3K4me3_Control.chr19.bam"
## download.file(H3K4me3_Control_BAM_url,"./H3K4me3_Control.chr19.bam")
## # get also the index file for previous file
## H3K4me3_Control_BAI_url="http://www.ebi.ac.uk/~ernesto/H3K4me3_Control.chr19.bam.bai"
## download.file(H3K4me3_Control_BAI_url,"./H3K4me3_Control.chr19.bam.bai")


###################################################
### code chunk number 25: example-cov.interval-histone-2 (eval = FALSE)
###################################################
## trm<-CoverageBamFile("./H3K4me3.chr19.bam",reads_mapped=864924)
## ctl<-CoverageBamFile("./H3K4me3_Control.chr19.bam",reads_mapped=319369)


###################################################
### code chunk number 26: example-cov.interval-histone-3 (eval = FALSE)
###################################################
## chr19_246050_247000.ratios=cov.interval(trm,ctl,bin_width=1,chr="chr19",
## start=246050,end=247000,do='ratio')


###################################################
### code chunk number 27: example-cov.interval-histone-4 (eval = FALSE)
###################################################
## export.wig(chr19_246050_247000.ratios,outfile="chr19_246050_247000.ratios.wig")


###################################################
### code chunk number 28: example-cov.interval-histone-5 (eval = FALSE)
###################################################
## chr19.ratios=cov.interval(trm,ctl,bin_width=1,chr="chr19",do='ratio')


