# Function to export the coverage information to a WIG file

export.wig <- function(cov.vector = "numeric", outfile = "character") {
    # get attributes from cov.vector
    chr <- attr(cov.vector, "chr")
    start <- attr(cov.vector, "start")
    bin_width <- attr(cov.vector, "bin_width")
    str <- c()
    if (bin_width > 1) {
        # track definition in order to cope with the binned data
        str <- paste("track type=wiggle_0 name=\"Rtrack\"\nfixedStep chrom=", chr, 
            " start=", start, " step=", bin_width, " span=", bin_width, sep = "")
    } else {
        str <- paste("track type=wiggle_0 name=\"Rtrack\"\nfixedStep chrom=", chr, 
            " start=", start, " step=1", sep = "")
    }
    write(str, outfile)
    write(cov.vector, outfile, ncolumns = 1, append = TRUE)
} 
