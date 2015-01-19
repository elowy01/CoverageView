# constructor
CoverageBamFile <- function(..., reads_mapped = 0L, run_type = c("single", "paired")) {
    reads_mapped <- as.integer(reads_mapped)
    run_type <- match.arg(run_type)
    .CoverageBamFile(BamFile(...), reads_mapped = reads_mapped, run_type = run_type)
}

setMethod(show, "CoverageBamFile", function(object) {
    callNextMethod()
    cat("reads_mapped:", object$reads_mapped, "\n")
    cat("run_type:", object$run_type, "\n")
}) 
