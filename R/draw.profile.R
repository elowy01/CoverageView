### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### method to draw a coverage profile

# draw profile for a single matrix
setMethod("draw.profile", signature("matrix"), function(data, outfile, ...) {
    # check arguments
    Lst <- list(...)
    if (is.null(Lst$ylab)) {
        stop("Error! Provide a value for 'ylab' parameter...")
    }
    # calculate mean coverage per row (position)
    v1 <- apply(data, 1, mean)
    v1 <- c(v1[1], v1)
    DF2 <- c()
    if (!is.null(attributes(data)$extend)) {
        extend <- attr(data, "extend")
        DF2 <- cbind(seq(from = (0 - extend), to = (0 + extend), by = extend * 2/(length(v1) - 
            1)), v1)
    } else if (!is.null(attributes(data)$no_windows)) {
        if (attributes(data)$offset > 0) {
            offset <- attributes(data)$offset
            by <- 100/attributes(data)$no_windows
            this_from <- -offset * by
            this_to <- 100 + offset * by
            DF2 <- cbind(seq(from = this_from, to = this_to, by = by), v1)
        } else {
            DF2 <- cbind(seq(from = 0, to = 100, by = 100/attributes(data)$no_windows), 
                v1)
        }
    }
    # if sys.parent=0 then the caller is not draw.heatmap
    if (sys.parent() == 0) {
        png(outfile, width = 700, height = 700)
    }
    # set default values for xlab
    if (is.null(Lst$xlab)) {
        if (!is.null(attributes(data)$extend)) {
            Lst$xlab <- "coord(bp)"
        } else if (!is.null(attributes(data)$no_windows)) {
            Lst$xlab <- "% gene"
        }
    }
    Lst$x <- DF2
    Lst$type <- "l"
    Lst$axes <- "F"
    # plot function
    do.call("plot", Lst)
    # add vertical lines at 0 an 100
    if (!is.null(attributes(data)$no_windows)) {
        if (attributes(data)$offset > 0) {
            abline(v = 0, col = 2, lty = 3, lwd = 2)
            abline(v = 100, col = 2, lty = 3, lwd = 2)
        }
    }
    grid()
    box()
    # set cex.laxis to 1.2 if cex is not passed as an argument
    if (is.null(Lst$cex.axis)) {
        Lst$cex.axis <- 1.2
    }
    par(cex.axis = Lst$cex.axis)
    # X-axis
    if (!is.null(attributes(data)$extend)) {
        # Calculate X tickmark interval
        interval <- (extend - (-extend))/10
        axis(1, at = seq(-extend, extend, interval))
    } else if (!is.null(attributes(data)$no_windows)) {
        if (attributes(data)$offset > 0) {
            axis(1, at = seq(0, 100, 10), cex.axis = Lst$cex.axis)
        } else {
            axis(1, at = seq(0, 100, 10), cex.axis = Lst$cex.axis)
        }
    }
    # Y-axis
    axis(2, las = 1)
    if (sys.parent() == 0) {
        dev.off()
    }
})

# draw profile for a list of data.frames
setMethod("draw.profile", signature("list"), function(data, outfile, ...) {
    # check arguments
    Lst <- list(...)
    if (is.null(Lst$ylab)) {
        stop("Error! Provide a value for 'ylab' parameter...")
    }
    # check if discrepancies among the attributes in the matrices
    no_windows_attr_v <- unlist(lapply(data, function(q) attributes(q)$no_windows))
    extend_attr_v <- unlist(lapply(data, function(q) attributes(q)$extend))
    if (length(no_windows_attr_v) > 1) {
        if (!isTRUE(all.equal(max(no_windows_attr_v), min(no_windows_attr_v)))) {
            stop("[ERROR]. no_windows argument value is different among the different coverage matrices, please rerun the cov.matrix function with the same no_windows value for all objects!")
        }
        offset_attr_v <- unlist(lapply(data, function(q) attributes(q)$offset))
        if (!isTRUE(all.equal(max(offset_attr_v), min(offset_attr_v)))) {
            stop("[ERROR]. offset argument value is different among the different coverage matrices, please rerun the cov.matrix function with the same offset value for all objects!")
        }
    } else if (length(extend_attr_v) > 1) {
        if (!isTRUE(all.equal(max(extend_attr_v), min(extend_attr_v)))) {
            stop("[ERROR]. extend argument value is different among the different coverage matrices, please rerun the cov.matrix function with the same extend value for all objects!")
        }
    }
    # if sys.parent=0 then the caller is not draw.heatmap
    if (sys.parent() == 0) {
        png(outfile, width = 700, height = 700)
    }
    legendnames <- c()
    # generate xlab depending on the extend/no_windows attributes
    this_xlab <- c()
    # max value for Y-axis
    maximum <- 0
    # iterate over each matrix
    data.post <- lapply(data, function(q) {
        # create legend names using filenames in 'data'
        legendnames <<- c(legendnames, attr(q, "filename"))
        v1 <- apply(q, 1, mean)
        v1 <- c(v1[1], v1)
        if (max(v1) > maximum) {
            maximum <<- c(max(v1))
        }
        DF2 <- c()
        if (!is.null(attributes(q)$extend)) {
            extend <- attr(q, "extend")
            DF2 <- cbind(seq(from = (0 - extend), to = (0 + extend), by = extend * 
                2/(length(v1) - 1)), v1)
            this_xlab <<- "coord(bp)"
        } else if (!is.null(attributes(q)$no_windows)) {
            if (attributes(q)$offset > 0) {
                offset <- attributes(q)$offset
                by <- 100/attributes(q)$no_windows
                this_from <- -offset * by
                this_to <- 100 + offset * by
                DF2 <- cbind(seq(from = this_from, to = this_to, by = by), v1)
            } else {
                DF2 <- cbind(seq(from = 0, to = 100, by = 100/attributes(q)$no_windows), 
                  v1)
            }
            this_xlab <<- "% gene"
        }
        return(DF2)
    })
    # plot the data.frames using a lapply
    color <- 0
    lapply(data.post, function(q) {
        color <<- color + 1
        plot(q, type = "l", axes = FALSE, col = color, cex.axis = 1.2, xlab = this_xlab, 
            ylim = c(0, round(maximum, digits = 1)), ...)
        par(new = TRUE)
    })
    # add vertical lines at 0 an 100
    if (!is.null(attributes(data[[1]])$no_windows)) {
        if (attributes(data[[1]])$offset > 0) {
            abline(v = 0, col = 2, lty = 3, lwd = 2)
            abline(v = 100, col = 2, lty = 3, lwd = 2)
        }
    }
    grid()
    box()
    # X-axis
    if (!is.null(attributes(data[[1]])$extend)) {
        # Calculate X tickmark interval
        extend <- attr(data[[1]], "extend")
        interval <- (extend - (-extend))/10
        axis(1, at = seq(-extend, extend, interval), cex.axis = 1)
    } else if (!is.null(attributes(data[[1]])$no_windows)) {
        if (attributes(data[[1]])$offset > 0) {
            axis(1, at = seq(0, 100, 10), cex.axis = 0.6)
        } else {
            axis(1, at = seq(0, 100, 10), cex.axis = 1)
        }
    }
    # Y-axis
    axis(2, las = 1, cex.axis = 1)
    # add legend
    legend("topright", legendnames, lty = 1, col = (1:color), cex = 1.2)
    if (sys.parent() == 0) {
        dev.off()
    }
}) 
