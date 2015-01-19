### method to draw a heatmap with the coverage levels

setMethod("draw.heatmap", signature("matrix"), function(data, outfile, color, ...) {
    # get the 0.90 quantile
    quantile <- quantile(as.vector(data), 0.9)
    # substitute all coverages greater than the 0.90 quantile by the coverage of the
    # 0.90 quantile in order to clean the outliers
    m1 <- apply(data, c(1, 2), function(x) if (x > quantile) 
        x <- quantile else x <- x)
    min <- min(m1)
    max <- max(m1)
    # calculate the break points for each color it works better than the first and
    # last % break is 0% and 100%
    this_breaks <- sapply(c(0, 65, 90, 100), function(x) max * (x)/100)
    # create a data.frame in which rows are the different genes and columns are the
    # different positions
    DFm1 <- as.data.frame(t(m1))
    # calculate the average coverage per row and add it as an additional column
    DFm1$avg <- apply(DFm1, 1, mean)
    # sort by avg using the reverse order
    DFm1_sorted <- DFm1[with(DFm1, (order(DFm1$avg))), ]
    # drop the $avg column
    DFm1_sorted_b <- subset(DFm1_sorted, select = -c(ncol(DFm1_sorted)))
    colorRamp <- colorRampPalette(c("white", color))(3)
    # plot heatmap
    png(outfile, width = 400, height = 1000)
    par(new = FALSE, mgp = c(0.5, 0.3, 0), fig = c(0, 0.9, 0.35, 1), las = 2, tcl = -0.3)
    # heatmap par
    par(mar = c(0.3, 5, 3, 0), cex.lab = 1.5)
    image(0:ncol(DFm1_sorted_b), 0:nrow(DFm1_sorted_b), t(DFm1_sorted_b), col = colorRamp, 
        axes = FALSE, breaks = this_breaks, ..., xlab = "", ylab = "genes", cex.main = 2)
    par(new = TRUE, fig = c(0, 0.9, 0.06, 0.35), mar = c(5, 5, 1, 0), mgp = c(3.5, 
        0.3, 0))
    # generate profile of the avg coverage per position pass attributes from original
    # object
    attributes(m1) <- attributes(data)
    draw.profile(m1, ylab = "avg coverage", cex.lab = 2, cex.axis = 1.6)
    # scale par
    par(new = TRUE, fig = c(0.01, 0.2, 0.03, 0.04), mar = c(0.18, 0.2, 0, 0))
    imageScale(DFm1_sorted_b, col = colorRamp, breaks = this_breaks, xlab = "", ylab = "", 
        cex.axis = 1.2)
    dev.off()
})

# draw heat.map for a list of data.frames
setMethod("draw.heatmap", signature("list"), function(data, outfile, color, ...) {
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
    # plot heatmap
    png(outfile, width = 400 * length(data), height = 700)
    # accommodate more than one cov histogram per page
    layout(matrix(seq(from = 1, to = length(data) * 3, by = 1), 3, length(data), 
        byrow = FALSE), heights = c(5, 1.5, 1.2))
    # layout(matrix(seq(from=1,to=length(data)*3,by=1), 3, length(data), byrow =
    # FALSE),heights=c(5,1,1.5))
    color <- 0
    # calculate the 0.90 quantile among all the DF within data
    quantile <- quantile(unlist(data), 0.9)
    data.post <- lapply(data, function(q) {
        # increment color value in order to assign a different color to each heatmap
        color <<- color + 1
        # substitute all coverages greater than the 0.90 quantile by the coverage of the
        # 0.90 quantile in order to clean the outliers
        m1 <- apply(q, c(1, 2), function(x) if (x > quantile) 
            x <- quantile else x <- x)
        min <- min(m1)
        max <- max(m1)
        # calculate the break points for each color it works better than the first and
        # last % break is 0% and 100%
        this_breaks <- sapply(c(0, 65, 90, 100), function(x) max * (x)/100)
        # create a data.frame in which rows are the different genes and columns are the
        # different positions
        DFm1 <- as.data.frame(t(m1))
        # calculate the average coverage per row and add it as an additional column
        DFm1$avg <- apply(DFm1, 1, mean)
        # sort by avg using the reverse order
        DFm1_sorted <- DFm1[with(DFm1, (order(DFm1$avg))), ]
        # drop the $avg column
        DFm1_sorted_b <- subset(DFm1_sorted, select = -c(ncol(DFm1_sorted)))
        colorRamp <- colorRampPalette(c("white", color))(3)
        par(new = FALSE, mgp = c(0, 0, 0), cex.axis = 0.5, las = 2, tcl = -0.3, mar = c(0, 
            5, 1, 1), cex.lab = 1.5)
        # heatmap par
        image(0:ncol(DFm1_sorted_b), 0:nrow(DFm1_sorted_b), t(DFm1_sorted_b), col = colorRamp, 
            axes = FALSE, breaks = this_breaks, ..., xlab = "", ylab = "genes", cex.main = 2)
        box()
        # create legend with filename
        legend("topleft", attr(q, "filename"), lty = 1, col = color, cex = 1.2)
        par(new = FALSE, mar = c(7, 5, 0, 1), mgp = c(4, 0.3, 0))
        # generate profile of the avg coverage per position pass attributes from original
        # object
        attributes(m1) <- attributes(q)
        draw.profile(m1, ylab = "avg coverage", cex.axis = 1.5)
        # scale par
        par(mar = c(8, 5, 2.6, 6))
        imageScale(DFm1_sorted_b, col = colorRamp, breaks = this_breaks, xlab = "", 
            ylab = "", cex.axis = 1.5)
    })
    dev.off()
}) 
