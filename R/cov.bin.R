### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### methods for calculating different types of coverage matrices with coverage per
### bin

# BAM files

# cov.bin relative to an anchor coord
setMethod("cov.bin", signature(extend = "numeric", no_windows = "missing", obj = "CoverageBamFile", 
    chr = "missing", grFromWig = "missing"), function(x, extend, obj, bin_width) {
    chr <- as.character(x[[1]])
    pos <- as.integer(x[[3]])
    strand <- as.character(x[[5]])
    message("[INFO] processing coords:", chr, " ", pos, " ", strand)
    # Create a GRange with the extended position
    rleToAdd <- Rle()
    min <- pos - extend
    max <- pos + (extend - 1)
    ceil <- max - min
    if (min <= 0) {
        # generate Rle object of the same length than the regions that is below 0
        rleToAdd <- Rle(rep(c(0), abs(min - 1)))
        # when min coordinate is below 1 then min will be equal to 1
        min <- 1
    }
    gr <- GRanges(seqnames = Rle(c(chr), c(1)), ranges = IRanges(start = min, end = max))
    param <- ScanBamParam(which = gr)
    # Read-in the BAM file
    reads <- c()
    # behave differently if type='paired'
    if (obj$run_type == "single") {
        reads <- readGAlignmentsFromBam(obj, param = param)
    } else {
        reads <- readGAlignmentPairsFromBam(obj, param = param)
    }
    # coverage method returns a Rle list
    cov <- coverage(reads)
    # get Rle for 'chr' only
    rleO <- cov[[chr]]
    # check if 'pos' is greater than the chromosome length
    if (pos > length(rleO)) {
        stop("[ERROR] processing coords:", chr, " ", pos, " ", strand, ".This position is greater than the length of the chromosome!")
    }
    rleSlice <- rleO[min:max]
    # if 'min' was below 0 then add a rle with 0 values just at the beggining of the
    # Rle object and a length equal to the number of nucleotides below 0 after
    # calculating the 'min' coordinate
    if (length(rleToAdd) > 1) {
        rleSlice <- append(rleToAdd, rleSlice)
    }
    # add a pseudocount of 1 to avoid divisions by 0
    rleSlice <- rleSlice + 1
    # initializing vectors for final rleSlice
    rleSlice_f <- c()
    # reverse rleSlice if '-' strand
    if (strand == "-") {
        rleSlice_f <- rev(rleSlice)
    } else {
        rleSlice_f <- rleSlice
    }
    if (bin_width == 1) {
        as.vector(rleSlice_f)
    } else {
        # calculate the break points for the aggregate function
        xstart <- seq(1, ceil, bin_width)
        xend <- seq(0 + bin_width, ceil + 1, bin_width)
        m <- aggregate(rleSlice_f, FUN = mean, start = xstart, end = xend)
        as.vector(m)
    }
})

# cov.bin for interval defined by chr,start,end
setMethod("cov.bin", signature(x = "missing", extend = "missing", no_windows = "missing", 
    obj = "CoverageBamFile", chr = "character", grFromWig = "missing"), function(obj, 
    chr, bin_width, start, end) {
    gr <- GRanges(seqnames = Rle(c(chr), c(1)), ranges = IRanges(start = start, end = end))
    param <- ScanBamParam(which = gr)
    # Read-in from the Treatment BAM file
    reads <- c()
    # behave differently if type='paired'
    if (obj$run_type == "single") {
        reads <- readGAlignmentsFromBam(obj, param = param)
    } else {
        reads <- readGAlignmentPairsFromBam(obj, param = param)
    }
    # coverage method returns a Rle list
    cov <- coverage(reads)
    # get Rle for 'this_chr' only
    rleO <- cov[[chr]]
    rleSlice <- rleO[start:end]
    if (bin_width == 1) {
        as.vector(rleSlice)
    } else {
        ceil <- end - start
        xstart <- seq(1, ceil, bin_width)
        xend <- seq(0 + bin_width, ceil + 1, bin_width)
        # correct incorrect split of the interval
        if (length(xstart) != length(xend)) {
            xstart <- xstart[1:length(xstart) - 1]
        }
        m <- aggregate(rleSlice, FUN = mean, start = xstart, end = xend)
        as.vector(m)
    }
})

# cov.bin for calculating the coverage from start to end, intervals defined in a
# file
setMethod("cov.bin", signature(extend = "missing", no_windows = "numeric", obj = "CoverageBamFile", 
    chr = "missing", grFromWig = "missing"), function(x, no_windows, obj, offset) {
    chr <- as.character(x[[1]])
    start.pre <- as.integer(x[[2]])
    end.pre <- as.integer(x[[3]])
    strand <- as.character(x[[5]])
    message("[INFO] processing coords:", chr, " ", start.pre, " ", end.pre)
    ceil <- end.pre - start.pre
    if (ceil < no_windows) {
        errormsg <- paste("[ERROR] nucleotidic interval from start to end coordinates cannot be smaller than the no_windows value. The interval:", 
            chr, start.pre, end.pre, "is not valid!", sep = " ")
        stop(errormsg, call. = FALSE)
    }
    # calculate the bin_width
    bin_width <- (ceil/no_windows)
    # extend interval if offset>0
    start.new <- 0
    end.new <- 0
    if (offset > 0) {
        start.new <- (start.pre - (offset * bin_width))
        end.new <- (end.pre + (offset * bin_width))
        ceil <- end.new - start.new
    } else {
        start.new <- start.pre
        end.new <- end.pre
    }
    if (start.new < 0) {
        stop("[ERROR] processing coords:", chr, " ", start.pre, " ", end.pre, ".After adjusting the start coordinate depending on the offset parameter value the new start coordinate is less than 0!. Please set the offset argument to 0")
    }
    # Create a GRange with the interval
    gr <- GRanges(seqnames = Rle(c(chr), c(1)), ranges = IRanges(start = start.new, 
        end = end.new))
    param <- ScanBamParam(which = gr)
    # Read-in from the Treatment BAM file
    reads <- c()
    # behave differently if type='paired'
    if (obj$run_type == "single") {
        reads <- readGAlignmentsFromBam(obj, param = param)
    } else {
        reads <- readGAlignmentPairsFromBam(obj, param = param)
    }
    # coverage method returns a Rle list
    cov <- coverage(reads)
    # get Rle for 'this_chr' only
    rleO <- cov[[chr]]
    # check if 'pos' is greater than the chromosome length
    if (start.new > length(rleO) || end.new > length(rleO)) {
        stop("[ERROR] processing coords:", chr, " ", start.pre, " ", end.pre, ".This position is greater than the length of the chromosome!. Try setting the offset parameter to 0")
    }
    rleSlice <- rleO[start.new:end.new]
    # initializing vectors for final rleSlice
    rleSlice_f <- c()
    # reverse rleSlice if '-' strand
    if (strand == "-") {
        rleSlice_f <- rev(rleSlice)
    } else {
        rleSlice_f <- rleSlice
    }
    # calculate the break points for the aggregate function
    xstart <- seq(from = 1, to = ceil, by = bin_width)
    xend <- seq(from = (0.9 + bin_width), to = ceil, by = bin_width)
    # check if lengths are correct and correct
    if (length(xend) != (no_windows + offset * 2)) {
        if (xend[length(xend)] <= ceil) {
            xend <- append(xend, xend[length(xend)] + bin_width)
        }
    }
    if (length(xstart) != (no_windows + offset * 2)) {
        if (xstart[length(xstart)] < ceil) {
            xstart <- append(xstart, xstart[length(xstart)] + bin_width)
        }
    }
    m <- aggregate(rleSlice_f, FUN = mean, start = xstart, end = xend)
    as.vector(m)
})

# WIG files

.cov.bin <- function(x, data, bin_width) {
    this_start <- x
    this_end <- (x + bin_width) - 1
    subrange <- window(data, start = this_start, end = this_end)
    return(mean(subrange))
}

setMethod("cov.bin", signature(extend = "numeric", no_windows = "missing", obj = "CoverageBigWigFile", 
    chr = "missing", grFromWig = "missing"), function(x, extend, obj, bin_width) {
    chr <- as.character(x[[1]])
    pos <- as.integer(x[[3]])
    strand <- as.character(x[[5]])
    message("[INFO] processing coords:", chr, " ", pos, " ", strand)
    min <- pos - extend
    max <- pos + extend
    # GRange with coordinates to get from WIG/BigWig file
    which <- GRanges(chr, IRanges(c(min), c(max - 1)))
    subGRange <- import(obj, which = which, asRangedData = FALSE,as="NumericList")
    covdata<-subGRange[[1]]
    
    # reverse if negative strand
    if (strand == "-") {
        covdata <- rev(covdata)
    }
    if (bin_width == 1) {
        # add a pseudocount of one to avoid divisions by 0
        covdata <- covdata + 1
        as.matrix(covdata)
    } else {
        # generate vector with the initial position of each bin
        start <- seq(1, extend * 2, bin_width)
        mean.cov.bin <- sapply(start, .cov.bin, data = covdata, bin = bin_width)
    }
})

setMethod("cov.bin", signature(extend = "missing", no_windows = "numeric", obj = "CoverageBigWigFile", 
    chr = "missing", grFromWig = "missing"), function(x, no_windows, obj, offset) {
    chr <- as.character(x[[1]])
    this_start <- as.integer(x[[2]])
    this_end <- as.integer(x[[3]])
    strand <- as.character(x[[5]])
    message("[INFO] processing coords:", chr, " ", this_start, " ", this_end)
    ceil <- this_end - this_start
    bin_width <- (ceil/no_windows)
    # extend interval if offset>0
    if (offset > 0) {
        if ((this_start - (offset * bin_width)) > 0) {
            this_start <- this_start - (offset * bin_width)
            this_end <- this_end + (offset * bin_width)
            ceil <- this_end - this_start
        } else {
            stop("[ERROR] processing coords:", chr, " ", this_start, " ", this_end, 
                ".After adjusting the start coordinate depending on the offset parameter value the new start coordinate is less than 0!. Please set the offset argument to 0")
        }
    }
    if (ceil < no_windows) {
        errormsg <- paste("[ERROR] nucleotidic interval from start to end coordinates cannot be smaller than the no_windows value. The interval:", 
            chr, this_start, this_end, "is not valid!", sep = " ")
        stop(errormsg)
    }
    
    # GRange with coordinates to get from WIG/BigWig file
    which <- GRanges(chr, IRanges(c(this_start), c(this_end)))
    subGRange <- import(obj, which = which, asRangedData = FALSE,as="NumericList")
    covdata<-subGRange[[1]]
    # reverse if negative strand
    if (strand == "-") {
        covdata <- rev(covdata)
    }
    # generate vector with the initial position of each bin
    start_list <- seq(1, ceil, bin_width)
    mean.cov.bin <- sapply(start_list, .cov.bin, data = covdata, bin = bin_width)
    if (offset > 0) {
        if ((this_start - (offset * bin_width)) < 0) {
            mean.cov.bin <- c(rep(mean(mean.cov.bin), offset), mean.cov.bin, rep(mean(mean.cov.bin), 
                offset))
        }
    }
    return(mean.cov.bin)
})

# cov.bin for interval defined by chr,start,end
setMethod("cov.bin", signature(x = "missing", extend = "missing", no_windows = "missing", 
    grFromWig = "GRanges", obj = "missing", chr = "character"), function(grFromWig, 
    chr, bin_width, start, end) {
    # split by chromosomes
    grFromWig_chr <- split(grFromWig, seqnames(grFromWig))
    # create sub GRanges from min to max
    subGRange <- grFromWig_chr[[chr]][start:(end - 1)]
    if (bin_width == 1) {
        sGR.DF <- subGRange@elementMetadata
        scores <- sGR.DF$score
        # add a pseudocount of one to avoid divisions by 0
        scores <- scores + 1
        as.matrix(scores)
    } else {
        ceil <- end - start
        start <- seq(1, ceil, bin_width)
        mean.cov.bin <- sapply(start, .cov.bin, data = subGRange, bin = bin_width)
    }
}) 
