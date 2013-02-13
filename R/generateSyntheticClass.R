generateSyntheticClass <- function(x, ...) {

    # Check arguments ----------------------------------------------------------
    
    # Check x.
    if (!(class(x) %in% c("big.matrix", "matrix", "data.frame"))) {
        stop("Argument x must be a big.matrix, matrix or data.frame.")
    }
    
    
    
    # Initialize ---------------------------------------------------------------
    xrows <- nrow(x)
    xcols <- ncol(x)
    xclass <- class(x)
    orgrows <- seq_len(xrows)
    newrows <- (xrows + 1L):(2L * xrows)
    
    
    
    # Set up x -----------------------------------------------------------------
    
    xold <- x
    if (xclass == "big.matrix") {
        x <- big.matrix(2 * xrows, xcols, type=typeof(xold),
                        dimnames=list(NULL, dimnames(x)[[2]]), ...)
        
        for (j in seq_len(xcols)) {
            x[orgrows, j] <- xold[orgrows, j]
        }
    } else if (xclass == "matrix") {
        x <- rbind(xold, matrix(vector(mode=typeof(x)), xrows, xcols))
    } else if (xclass == "data.frame") {
        x <- rbind(xold, xold)
    }
    
    
    
    # Synthesize new rows and y ------------------------------------------------
    
    for (j in seq_len(xcols)) {
        x[newrows, j] <- xold[as.integer(runif(xrows, 0, xrows) + 1), j]
    }
    
    y <- factor(c(rep.int(1L, xrows), rep.int(2L, xrows)),
                labels=c("Original", "Synthetic"))
    
    return(list(x=x, y=y))
}
