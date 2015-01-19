setGeneric("cov.bin", function(x, extend, no_windows, grFromWig, obj, chr, bin_width, 
    offset, start, end) standardGeneric("cov.bin"))
setGeneric("cov.interval", function(tr, ctl, normalization = NULL, chr, start, end, 
    bin_width = 1, do = "log2ratio") standardGeneric("cov.interval"))
setGeneric("cov.matrix", function(tr, ctl, coordfile, normalization = NULL, bin_width = 1, 
    extend, no_windows, offset = 0, do = "log2ratio", num_cores) standardGeneric("cov.matrix"))
setGeneric("draw.profile", function(data, outfile = NULL, ...) standardGeneric("draw.profile"))
setGeneric("draw.heatmap", function(data, outfile = NULL, color = "blue", ...) standardGeneric("draw.heatmap"))
setGeneric("genome.covplot.depth", function(data, outfile = NULL, max_depth = 50) standardGeneric("genome.covplot.depth"))
setGeneric("genome.covplot.cumdepth", function(data, outfile = NULL, max_depth = 20) standardGeneric("genome.covplot.cumdepth"))
setGeneric("write.profile", function(DF, outfile = NULL) standardGeneric("write.profile")) 
