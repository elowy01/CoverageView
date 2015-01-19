### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### methods for generating a coverage profile for a certain genomic interval

### cov.interval for two BAM files
setMethod("cov.interval", signature(tr = "CoverageBamFile", ctl = "CoverageBamFile"), 
    function(tr, ctl, normalization, chr, start, end, bin_width, do) {
        # check validity of the 'do' argument
        if (do != "log2ratio" && do != "ratio" && do != "substraction") {
            stop("[ERROR] 'do' argument value is not correct. Possible values are 'log2ratio'/'ratio'/'substraction'")
        }
        this_start <- c()
        this_end <- c()
        if (!missing(start) && !missing(end)) {
            # check validity of the 'start' and 'end' args
            if (start < 1 || end < start) {
                stop("[ERROR] Either the start or end argument is not valid, please check them!")
            }
            # round to the nearest ten
            this_start <- round(start - 4, -1)
            if (this_start == 0) {
                this_start <- 1
            }
            this_end <- end
        } else {
            # if 'chr' is set and neither the 'start' and 'end' args, then the operation will
            # be performed over the entire chromosome
            this_start <- 1
            this_end <- scanBamHeader(tr$path)[[1]][["targets"]][chr]
        }
        this_chr <- chr
        treat.vector <- cov.bin(obj = tr, chr = chr, start = this_start, end = this_end, 
            bin_width = bin_width)
        ctrl.vector <- cov.bin(chr = chr, start = this_start, end = this_end, obj = ctl, 
            bin_width = bin_width)
        if (!is.null(normalization)) {
            if (normalization == "rpm") {
                if (tr$reads_mapped == 0) {
                  message("[INFO] No readcount provided for ", deparse(substitute(tr)), 
                    " object.Starting readcount...")
                  # count number of mapped reads in the BAM file (for normalization)
                  param <- ScanBamParam()
                  bamFlag(param) <- scanBamFlag(isUnmappedQuery = FALSE, isNotPrimaryRead = FALSE)
                  c <- countBam(tr, param = param)
                  readcount <- c$records
                  msg1 <- paste("[INFO] Number of mapped reads is", readcount, "...")
                  message(msg1)
                  tr$reads_mapped <- readcount
                }
                # normalize
                treat.vector <- .do.normalization(treat.vector, tr$reads_mapped)
                if (ctl$reads_mapped == 0) {
                  message("[INFO] No readcount provided for ", deparse(substitute(ctl)), 
                    " object.Starting readcount...")
                  # count number of mapped reads in the BAM file (for normalization)
                  bamFlag(param) <- scanBamFlag(isUnmappedQuery = FALSE, isNotPrimaryRead = FALSE)
                  c <- countBam(ctl, param = param)
                  readcount <- c$records
                  msg1 <- paste("[INFO] Number of mapped reads is", readcount, "...")
                  message(msg1)
                  ctl$reads_mapped <- readcount
                }
                # normalize
                ctrl.vector <- .do.normalization(ctrl.vector, ctl$reads_mapped)
            }
        }
        computed.vector <- c()
        # add a pseudocount of 1 to avoid divisions by 0 if 'do' involves a ratio
        if (do == "log2ratio") {
            treat.vector <- treat.vector + 1
            ctrl.vector <- ctrl.vector + 1
            computed.vector <- log2(treat.vector/ctrl.vector)
        } else if (do == "ratio") {
            treat.vector <- treat.vector + 1
            ctrl.vector <- ctrl.vector + 1
            computed.vector <- (treat.vector/ctrl.vector)
        } else if (do == "substraction") {
            computed.vector <- (treat.vector - ctrl.vector)
        }
        # add attributes necessary for function export.wig
        attr(computed.vector, "start") <- this_start
        attr(computed.vector, "chr") <- chr
        attr(computed.vector, "bin_width") <- bin_width
        return(computed.vector)
    })

## cov.interval for two WIG files

setMethod("cov.interval", signature(tr = "CoverageBigWigFile", ctl = "CoverageBigWigFile"), 
    function(tr, ctl, normalization, chr, start, end, bin_width = 1, do) {
        # check validity of the 'do' argument
        if (do != "log2ratio" && do != "ratio" && do != "substraction") {
            stop("[ERROR] 'do' argument value is not correct. Possible values are 'log2ratio'/'ratio'/'substraction'")
        }
        this_chr <- chr
        # import data from wig/BigWig file
        treat.grFromWig <- import(tr, asRangedData = FALSE)
        treat.grFromWig <- as(treat.grFromWig, "GenomicRanges")
        con.grFromWig <- import(ctl, asRangedData = FALSE)
        con.grFromWig <- as(con.grFromWig, "GenomicRanges")
        ######## check the coordinates for which the profile will be generated
        this_start <- c()
        this_end <- c()
        if (!missing(start) && !missing(end)) {
            # check validity of the 'start' and 'end' args
            if (start < 1 || end < start) {
                stop("[ERROR] Either the start or end argument is not valid, please check them!")
            }
            # round to the nearest ten
            this_start <- round(start - 4, -1)
            if (this_start == 0) {
                this_start <- 1
            }
            this_end <- end
        } else {
            # if 'chr' is set and neither the 'start' and 'end' args, then the operation will
            # be performed over the entire chromosome
            this_start <- 1
            # estimate 'this_end' by calculating the length of the 'wig' file for this chr
            f <- seqnames(treat.grFromWig)
            this_end <- runLength(f[f == chr])
        }
        
        treat.vector <- cov.bin(chr = chr, start = this_start, end = this_end, bin_width = bin_width, 
            grFromWig = treat.grFromWig)
        ctrl.vector <- cov.bin(chr = chr, start = this_start, end = this_end, bin_width = bin_width, 
            grFromWig = con.grFromWig)
        if (!is.null(normalization)) {
            if (normalization == "rpm") {
                if (tr@reads_mapped == 0) {
                  stop("[ERROR] No readcount provided for ", deparse(substitute(ifile)), 
                    " object.Please provide this number...")
                }
                # normalize
                treat.vector <- .do.normalization(treat.vector, tr@reads_mapped)
                if (ctl@reads_mapped == 0) {
                  stop("[ERROR] No readcount provided for ", deparse(substitute(ifile)), 
                    " object.Please provide this number...")
                }
                # normalize
                ctrl.vector <- .do.normalization(ctrl.vector, ctl@reads_mapped)
            }
        }
        computed.vector <- c()
        # add a pseudocount of 1 to avoid divisions by 0 if 'do' involves a ratio
        if (do == "log2ratio") {
            treat.vector <- treat.vector + 1
            ctrl.vector <- ctrl.vector + 1
            computed.vector <- log2(treat.vector/ctrl.vector)
        } else if (do == "ratio") {
            treat.vector <- treat.vector + 1
            ctrl.vector <- ctrl.vector + 1
            computed.vector <- (treat.vector/ctrl.vector)
        } else if (do == "substraction") {
            computed.vector <- (treat.vector - ctrl.vector)
        }
        # add attributes necessary for function export.wig
        attr(computed.vector, "start") <- this_start
        attr(computed.vector, "chr") <- this_chr
        attr(computed.vector, "bin_width") <- bin_width
        return(computed.vector)
    }) 
