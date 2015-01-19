.do.normalization <- function(m, readcount) {
    # rpm normalization
    m <- (m/readcount) * 1e+06
    return(m)
} 
