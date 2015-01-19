### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### methods for generating different types of genome coverage plots

# Function to generate a plot showing the number of genomic positions reaching a
# certain read depth (for a single CoverageBamFile)
setMethod("genome.covplot.depth", signature(data = "CoverageBamFile"), function(data, 
    outfile, max_depth) {
    # Read-in the BAM file
    reads <- c()
    if (data$run_type == "single") {
        reads <- readGAlignmentsFromBam(data)
    } else {
        reads <- readGAlignmentPairsFromBam(data)
    }
    # coverage method returns a Rle list
    cov <- coverage(reads)
    tableSum <- lapply(cov, function(q) {
        t <- table(q)
        # consider only coverage read depths up to 50X
        subt <- t[as.character(0:max_depth)]
        v <- as.vector(subt)
        # replacing all NAs into 0s
        v[is.na(v)] <- 0
        return(v)
    })
    # calculate the sums for all the chromosomes
    tableSum2 <- Reduce("+", tableSum)
    png(outfile, width = 700, height = 700)
    par(mar = c(9, 9, 2, 2), mgp = c(7, 1, 0))
    plot(0:max_depth, tableSum2/1000, axes = FALSE, type = "b", cex.lab = 2, lty = 1, 
        xlab = "Read depth", ylab = "Number of bases (M)", col = "red")
    axis(1, cex.axis = 1.2)
    axis(2, cex.axis = 1.2, las = 2)
    dev.off()
})

# Function to estimate the number of positions with a certain read depth (for a
# list of CoverageBamFile objects)
setMethod("genome.covplot.depth", signature(data = "list"), function(data, outfile, 
    max_depth) {
    # iterate over initialize graphic paramaters for setting them iteratively
    png(outfile, width = 1000, height = 650)
    par(mar = c(9, 9, 2, 23), mgp = c(7, 1, 0), xpd = TRUE)
    color <- 0
    this_lty <- 0
    this_pch <- 0
    legendnames <- c()
    data.max <- lapply(data, function(q) {
        color <<- color + 1
        this_lty <<- this_lty + 1
        this_pch <<- this_pch + 1
        legendnames <<- c(legendnames, basename(q$path))
        # Read-in the BAM file
        reads <- c()
        if (q$run_type == "single") {
            reads <- readGAlignmentsFromBam(q)
        } else {
            reads <- readGAlignmentPairsFromBam(q)
        }
        # coverage method returns a Rle list
        cov <- coverage(reads)
        tableSum <- lapply(cov, function(q) {
            t <- table(q)
            # consider only coverage read depths up to 50X
            subt <- t[as.character(0:max_depth)]
            v <- as.vector(subt)
            # replacing all NAs into 0s
            v[is.na(v)] <- 0
            return(v)
        })
        # calculate the sums for all the chromosomes
        tableSum2 <- Reduce("+", tableSum)
        plot(0:max_depth, tableSum2/1000, axes = FALSE, type = "b", cex = 0.6, cex.lab = 1.5, 
            lty = this_lty, pch = this_pch, xlab = "Read depth", ylab = "Number of bases (M)", 
            col = color)
        par(new = TRUE)
    })
    axis(1, cex.axis = 1.5)
    axis(2, cex.axis = 1.5, las = 2)
    # add legend with negative inset in order to draw the legend out of the plotting
    # area
    legend("topright", inset = c(-0.6, 0), legendnames, lty = (1:length(data)), col = (1:length(data)), 
        pch = (1:length(data)), cex = 1.2)
    dev.off()
})

# Function to calculate the cumulative read depth (for a single CoverageBamFile
# object)
setMethod("genome.covplot.cumdepth", signature(data = "CoverageBamFile"), function(data, 
    outfile, max_depth) {
    # Read-in the BAM file
    reads <- c()
    if (data$run_type == "single") {
        reads <- readGAlignmentsFromBam(data)
    } else {
        reads <- readGAlignmentPairsFromBam(data)
    }
    # coverage method returns a Rle list
    cov <- coverage(reads)
    # max coverage value
    max_cov <- (max((max(cov)))) + 1
    # total genome size
    genome_size <- 0
    tableSum <- lapply(cov, function(q) {
        genome_size <<- genome_size + length(q)
        t <- table(q)
        # consider only coverage read depths up to the max coverage
        subt <- t[as.character(0:max_cov)]
        v <- as.vector(subt)
        # replacing all NAs into 0s
        v[is.na(v)] <- 0
        return(v)
    })
    # calculate the sums for all the chromosomes
    tableSum2 <- Reduce("+", tableSum)
    # calculate vector with cumulative read depth for depths from 0X to max_depth
    v <- sapply(1:(max_depth + 1), function(x) sum(tableSum2[x:max_cov]) * 100/genome_size)
    png(outfile, width = 700, height = 700)
    plot(c(0:max_depth), v, type = "b", xlim = c(0, max_depth), ylim = c(0, 100), 
        cex.lab = 1.5, lty = 1, las = 2, xlab = "Cumulative read depth", ylab = "Percentage of genome (%)", 
        col = "red", axes = FALSE)
    # draw axes
    par(mar = c(5, 5.5, 2, 2))
    axis(1, cex.axis = 1.2)
    axis(2, at = seq(from = 0, to = 100, by = 5), cex.axis = 1.2, las = 2)
    dev.off()
})

# Function to calculate the cumulative read depth (for a list of CoverageBamFile
# objects)
setMethod("genome.covplot.cumdepth", signature(data = "list"), function(data, outfile, 
    max_depth) {
    # iterate over initialize graphic paramaters for setting them iteratively
    png(outfile, width = 900, height = 650)
    # Add extra space to right of plot area; change clipping to figure
    par(mar = c(4.5, 4.5, 6, 22), xpd = TRUE)
    color <- 0
    this_lty <- 0
    this_pch <- 0
    legendnames <- c()
    data.max <- lapply(data, function(q) {
        color <<- color + 1
        this_lty <<- this_lty + 1
        this_pch <<- this_pch + 1
        legendnames <<- c(legendnames, basename(q$path))
        # Read-in the BAM file
        reads <- c()
        if (q$run_type == "single") {
            reads <- readGAlignmentsFromBam(q)
        } else {
            reads <- readGAlignmentPairsFromBam(q)
        }
        # coverage method returns a Rle list
        cov <- coverage(reads)
        # max coverage value
        max_cov <- (max((max(cov)))) + 1
        # calculate total genome size
        genome_size <- 0
        tableSum <- lapply(cov, function(q) {
            genome_size <<- genome_size + length(q)
            t <- table(q)
            # consider only coverage read depths up to the max coverage
            subt <- t[as.character(0:max_cov)]
            v <- as.vector(subt)
            # replacing all NAs into 0s
            v[is.na(v)] <- 0
            return(v)
        })
        # calculate the sums for all the chromosomes
        tableSum2 <- Reduce("+", tableSum)
        # calculate vector with cumulative read depth for depths from 1X to max_depth
        v <- sapply(1:(max_depth + 1), function(x) sum(tableSum2[x:max_cov]) * 100/genome_size)
        plot(c(0:max_depth), v, type = "b", ylim = c(0, 100), axes = FALSE, cex.lab = 1.5, 
            cex = 0.6, lty = this_lty, pch = this_pch, col = color, las = 2, xlab = "Cumulative read depth", 
            ylab = "Percentage of genome (%)")
        par(new = TRUE)
        return(range(v))
    })
    # X-axis
    axis(1, cex.axis = 1.2, las = 1)
    # Y-axis
    axis(2, at = seq(from = 0, to = 100, by = 5), cex.axis = 1.2, las = 2)
    # add legend out of the plotting area using negative inset
    legend("topright", inset = c(-0.6, 0), legendnames, lty = (1:length(data)), col = (1:length(data)), 
        pch = (1:length(data)), cex = 1.2)
    dev.off()
}) 
