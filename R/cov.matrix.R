### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### methods for generating different types of coverage profiles

### function to handle the errors propagated from cov.bin function
.checkErrors <- function(x) {
    if (class(x) == "try-error") {
        # die if any error
        stop(x[1])
    }
}

# helper function to set the num_cores option and set num_cores=1 if platform is
# Windows
.defaultCores <- function(x) {
    if (.Platform$OS.type == "windows") {
        return(1)
    } else {
        if (!missing(x)) {
            return(x)
        } else {
            return(detectCores())
        }
    }
}

### cov.matrix for a single BAM/BigWIG file

setMethod("cov.matrix", signature(tr = "ANY", ctl = "missing"), function(tr, coordfile, 
    normalization, bin_width, extend, no_windows, offset, num_cores) {
    message("[INFO] A coverage matrix composed of average coverages will be generated.")
    # get number of cores to use respecting the user's global setting
    num_cores <- getOption("mc.cores", .defaultCores(x = num_cores))
    # check validity of arguments
    if (missing(extend) && missing(no_windows)) {
        stop("[ERROR] Specify the 'extend' or 'no_windows' argument!")
    } else if (!missing(extend) && !missing(no_windows)) {
        stop("[ERROR] Both the 'extend' and 'no_windows' arguments cannot be set simultaneously, please set only one of both!")
    }
    if (bin_width < 1) {
        stop("[ERROR] bin_width parameter is not valid. Enter a different value")
    }
    # get extension of the file containing the coordinates from its path
    ext <- file_ext(coordfile)
    # create a GRanges object from the coordinate file
    grDF <- c()
    # check the format of the file and act accordingly
    if (length(grep("bed", ext, ignore.case = TRUE)) > 0) {
        # create a GRanges from a bed file
        grDF <- .bed.to.GRanges(coordfile)
    } else if (length(grep("gff", ext, ignore.case = TRUE)) > 0) {
        # create a GRanges from a gff file
        grDF <- .gff.to.GRanges(coordfile)
    }
    ## coverage relative to an anchor coord
    mres <- c()
    if (!missing(extend)) {
        # trick in order to use the mclapply
        loopfun <- function(i) {
            cov.bin(grDF[i, ], extend = extend, obj = tr, bin_width = bin_width)
        }
        # returns a list
        res <- mclapply(1:nrow(grDF), loopfun, mc.cores = num_cores)
        lapply(res, .checkErrors)
        # returns a matrix
        mres <- sapply(res, function(x) x)
    } else if (!missing(no_windows)) {
        ## coverage for intervals defined by a start/end coord trick in order to use the
        ## mlcapply
        if (bin_width > 1) {
            message(paste("[INFO] no_windows parameter was set. bin_width=", bin_width, 
                " parameter will be ignored in all the downstream analyses..."), 
                sep = "")
        }
        loopfun <- function(i) {
            cov.bin(grDF[i, ], obj = tr, no_windows = no_windows, offset = offset)
        }
        # returns a list
        res <- mclapply(1:nrow(grDF), loopfun, mc.cores = num_cores)
        # check if errors
        lapply(res, .checkErrors)
        # returns a matrix
        mres <- sapply(res, function(x) x)
    }
    if (!is.null(normalization)) {
        readcount <- 0
        if (normalization == "rpm") {
            if (class(tr) == "CoverageBamFile") {
                if (tr$reads_mapped == 0) {
                  msg <- paste("[INFO] No readcount provided for", deparse(substitute(tr)), 
                    "object.Starting readcount...")
                  message(msg)
                  # count number of mapped reads in the BAM file (for normalization)
                  param <- ScanBamParam()
                  bamFlag(param) <- scanBamFlag(isUnmappedQuery = FALSE, isNotPrimaryRead = FALSE)
                  c <- countBam(tr, param = param)
                  readcount <- c$records
                  msg1 <- paste("[INFO] Number of mapped reads is", readcount, "...")
                  message(msg1)
                } else {
                  readcount <- tr$reads_mapped
                }
            } else if (class(tr) == "CoverageBigWigFile") {
                if (tr@reads_mapped == 0) {
                  stop("[ERROR]. The 'reads_mapped' slot is not set for the CoverageBigWigFile object. \n                       The number of reads mapped is necessary for doing the 'rpm' normalization")
                } else {
                  readcount <- tr@reads_mapped
                }
            }
            mres <- .do.normalization(mres, readcount)
        }
    }
    # add an attribute with the filename
    attr(mres, "filename") <- basename(path(tr))
    if (!missing(extend)) {
        # add the 'extend' info into the data.frame
        attr(mres, "extend") <- extend
    } else if (!missing(no_windows)) {
        # add the 'no_windows' and 'offset' attributes into the data.frame
        attr(mres, "no_windows") <- no_windows
        attr(mres, "offset") <- offset
    }
    # returns a dataframe where the rows are the different positions/bins analyzed
    # and the columns are the different genes
    return(mres)
})

### cov.matrix for two BAM/WIG-BigWIG files that will be compared one against the
### other

setMethod("cov.matrix", signature(tr = "ANY", ctl = "ANY"), function(tr, ctl, coordfile, 
    normalization, bin_width, extend, no_windows, offset, do, num_cores) {
    # get number of cores to use respecting the user's global setting
    num_cores <- getOption("mc.cores", .defaultCores(x = num_cores))
    # check validity of arguments
    if (missing(extend) && missing(no_windows)) {
        stop("[ERROR] Specify the 'extend' or 'no_windows' argument!")
    } else if (!missing(extend) && !missing(no_windows)) {
        stop("[ERROR] Both the 'extend' and 'no_windows' arguments cannot be set simultaneously, please set only one of both!")
    }
    if (bin_width < 1) {
        stop("[ERROR] bin_width parameter is not valid. Enter a different value")
    }
    # check validity of the 'do' argument
    if (do != "log2ratio" && do != "ratio" && do != "substraction") {
        stop("[ERROR] 'do' argument value is not correct. Possible values are 'log2ratio'/'ratio'/'substraction'")
    }
    warningmsg <- paste("[INFO] A coverage matrix will be generated by calculating the '", 
        do, "' of the coverages.", sep = "")
    message(warningmsg)
    # get extension of the file containing the coordinates from its path
    ext <- file_ext(coordfile)
    # create a GRanges object from the coordinate file
    grDF <- c()
    # check the format of the file and act accordingly
    if (length(grep("bed", ext, ignore.case = TRUE)) > 0) {
        # create a GRanges from a bed file
        grDF <- .bed.to.GRanges(coordfile)
    } else if (length(grep("gff", ext, ignore.case = TRUE)) > 0) {
        # create a GRanges from a gff file
        grDF <- .gff.to.GRanges(coordfile)
    }
    mtreat_mat <- c()
    mctrl_mat <- c()
    if (!missing(extend)) {
        ## coverage relative to an anchor coord trick in order to use the mlcapply
        loopfun <- function(i, obj) {
            cov.bin(grDF[i, ], extend = extend, obj = obj, bin_width = bin_width)
        }
        
        ## treatment
        message("[INFO] processing:", path(tr))
        # returns a list
        treat_mat <- mclapply(1:nrow(grDF), loopfun, mc.cores = num_cores, obj = tr)
        # returns a matrix
        mtreat_mat <- sapply(treat_mat, function(x) x)
        
        ## control
        message("[INFO] processing:", path(ctl))
        # returns a list
        ctrl_mat <- mclapply(1:nrow(grDF), loopfun, mc.cores = num_cores, obj = ctl)
        # returns a matrix
        mctrl_mat <- sapply(ctrl_mat, function(x) x)
    } else if (!missing(no_windows)) {
        ## coverage for intervals defined by a start/end coord
        if (bin_width > 1) {
            message(paste("[INFO] no_windows parameter was set. bin_width=", bin_width, 
                " parameter will be ignored in all the downstream analyses..."), 
                sep = "")
        }
        # trick in order to use the mlcapply
        loopfun <- function(i, obj) {
            cov.bin(grDF[i, ], obj = tr, no_windows = no_windows, offset = offset)
            res <- cov.bin(grDF[i, ], obj = obj, bin_width = bin_width, no_windows = no_windows, 
                offset = offset)
            return(res)
        }
        
        ## treatment
        message("[INFO] processing:", path(tr))
        # returns a list
        treat_mat <- mclapply(1:nrow(grDF), loopfun, mc.cores = num_cores, obj = tr)
        # check if errors
        lapply(treat_mat, .checkErrors)
        # returns a matrix
        mtreat_mat <- sapply(treat_mat, function(x) x)
        
        ## control
        message("[INFO] processing:", path(ctl))
        # returns a list
        ctrl_mat <- mclapply(1:nrow(grDF), loopfun, mc.cores = num_cores, obj = ctl)
        # check if errors
        lapply(ctrl_mat, .checkErrors)
        # returns a matrix
        mctrl_mat <- sapply(ctrl_mat, function(x) x)
    }
    
    # normalization
    if (!is.null(normalization)) {
        if (normalization == "rpm") {
            readcount.trm <- 0
            readcount.ctl <- 0
            if (class(tr) == "CoverageBamFile") {
                if (tr$reads_mapped == 0) {
                  msg <- paste("No readcount provided for", deparse(substitute(tr)), 
                    "object.Starting readcount...")
                  message(msg)
                  # count number of mapped reads in the BAM file (for normalization)
                  param <- ScanBamParam()
                  bamFlag(param) <- scanBamFlag(isUnmappedQuery = FALSE, isNotPrimaryRead = FALSE)
                  c <- countBam(tr, param = param)
                  readcount.trm <- c$records
                  msg1 <- paste("[INFO] Number of mapped reads is", readcount.trm, 
                    "...")
                  message(msg1)
                }
            } else if (class(tr) == "CoverageBigWigFile") {
                if (tr@reads_mapped == 0) {
                  stop("[ERROR]. The 'reads_mapped' slot is not set for the CoverageBigWigFile object. \n                       The number of reads mapped is necessary for doing the 'rpm' normalization")
                } else {
                  readcount.trm <- tr@reads_mapped
                }
            }
            mtreat_mat <- .do.normalization(mtreat_mat, readcount.trm)
            if (class(ctl) == "CoverageBamFile") {
                if (ctl$reads_mapped == 0) {
                  msg <- paste("No readcount provided for", deparse(substitute(ctl)), 
                    "object.Starting readcount...")
                  message(msg)
                  # count number of mapped reads in the BAM file (for normalization)
                  param <- ScanBamParam()
                  bamFlag(param) <- scanBamFlag(isUnmappedQuery = FALSE, isNotPrimaryRead = FALSE)
                  c <- countBam(ctl, param = param)
                  readcount.ctl <- c$records
                  msg1 <- paste("[INFO] Number of mapped reads is", readcount.ctl, 
                    "...")
                  message(msg1)
                }
            } else if (class(ctl) == "CoverageBigWigFile") {
                if (ctl@reads_mapped == 0) {
                  stop("[ERROR]. The 'reads_mapped' slot is not set for the CoverageBigWigFile object. The number\n                    of reads mapped is necessary for doing the 'rpm' normalization")
                } else {
                  readcount.ctl <- ctl@reads_mapped
                }
            }
            mctrl_mat <- .do.normalization(mctrl_mat, readcount.ctl)
        }
    }
    
    operationDF <- c()
    # change operation to be performed depending on the 'do' value add a pseudocount
    # of 1 to avoid divisions by 0 if 'do' involves a ratio
    if (do == "log2ratio") {
        mtreat_mat <- mtreat_mat + 1
        mctrl_mat <- mctrl_mat + 1
        operationDF <- log2(mtreat_mat/mctrl_mat)
    } else if (do == "ratio") {
        mtreat_mat <- mtreat_mat + 1
        mctrl_mat <- mctrl_mat + 1
        operationDF <- (mtreat_mat/mctrl_mat)
    } else if (do == "substraction") {
        operationDF <- (mtreat_mat - mctrl_mat)
    }
    
    # add an attribute with the filename
    if (!missing(extend)) {
        # add the 'extend' info into the data.frame
        attr(operationDF, "extend") <- extend
    } else if (!missing(no_windows)) {
        # add the 'no_windows' attribute into the data.frame
        attr(operationDF, "no_windows") <- no_windows
        attr(operationDF, "offset") <- offset
    }
    return(operationDF)
}) 

