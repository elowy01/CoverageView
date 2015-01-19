# Allclasses

# CoverageBamFile inherits from Rsamtools BamFile class
.CoverageBamFile <- setRefClass("CoverageBamFile", fields = c(reads_mapped = "integer", 
    run_type = "character"), contains = "BamFile")


# CoverageBigWigFile inherits from rtracklayer BigWigFile class
.CoverageBigWigFile <- setClass("CoverageBigWigFile", representation = representation(reads_mapped = "integer"),
    contains = "BigWigFile") 
