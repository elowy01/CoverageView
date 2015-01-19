# method to write the coverage information in a file

setMethod("write.profile", signature("matrix"), function(DF, outfile) {
    # check arguments
    if (is.null(outfile)) {
        stop("Error! Provide a value for 'outfile' parameter...")
    }
    # calculate mean coverage per row (position)
    v1 <- apply(DF, 1, mean)
    v1 <- c(v1[1], v1)
    DF2 <- c()
    if (!is.null(attributes(DF)$extend)) {
        extend <- attr(DF, "extend")
        DF2 <- cbind(seq(from = (0 - extend), to = (0 + extend), by = extend * 2/(length(v1) - 
            1)), v1)
        colnames(DF2) <- c("coord", "cov")
    } else if (!is.null(attributes(DF)$no_windows)) {
        if (attributes(DF)$offset > 0) {
            offset <- attributes(DF)$offset
            by <- 100/attributes(DF)$no_windows
            this_from <- -offset * by
            this_to <- 100 + offset * by
            DF2 <- cbind(seq(from = this_from, to = this_to, by = by), v1)
        } else {
            DF2 <- cbind(seq(from = 0, to = 100, by = 100/attributes(DF)$no_windows), 
                v1)
        }
        colnames(DF2) <- c("%", "cov")
    }
    write.table(DF2, file = outfile, row.names = FALSE)
}) 
