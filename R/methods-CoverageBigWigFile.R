# constructor
CoverageBigWigFile <- function(..., reads_mapped = 0L) {
    reads_mapped <- as.integer(reads_mapped)
    .CoverageBigWigFile(BigWigFile(...), reads_mapped = reads_mapped)
}

setMethod(show, "CoverageBigWigFile", function(object) {
    callNextMethod()
    cat("reads_mapped:", object@reads_mapped, "\n")
}) 
