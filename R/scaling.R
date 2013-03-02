setGeneric("scaling", function(prox, ...) standardGeneric("scaling"))



setMethod("scaling", signature(prox="bigrfprox"),
          function(
    prox,
    nscale=2L,
    trace=0L) {
    
    # Check arguments ----------------------------------------------------------
    
    # Check trace.
    if (!is.numeric(trace) ||
            abs(trace - round(trace)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument trace must be an integer.")
    }
    trace <- as.integer(round(trace))
    if (trace < 0L || trace > 2L) {
        stop("Argument trace must be 0, 1 or 2.")
    }
    
    if (trace >= 1L) message("Checking arguments.")
    
    # Check prox.
    if (!class(prox) == "bigrfprox") {
        stop("Argument prox must be an object of class bigrfprox")
    }
    
    # Check nscale.
    if (!is.numeric(nscale) ||
            abs(nscale - round(nscale)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument nscale must be an integer.")
    }
    nscale <- as.integer(round(nscale))
    if (nscale < 1L) {
        stop("Argument trace must be at least 1.")
    }
    
    
    
    # Initialize ---------------------------------------------------------------
    
    if (trace >= 1L) message("Initializing.")
    
    nexamples <- nrow(prox)
    nnearest <- ncol(prox)
    scale <- matrix(numeric(), nexamples, nscale,
                    dimnames=list(Example=NULL, CoordNum=NULL))
    ev <- matrix(numeric(), nexamples, nscale)
    bl <- numeric(nscale)
    dl <- numeric(nscale)
    red <- foreach(i=seq_len(nexamples), .combine=c) %dopar% {
        sum(prox[i, ]) / nexamples
    }
    sred <- sum(red) / nexamples
    
    
    
    # Compute scaling coordinates ----------------------------------------------
    
    for (it in seq_len(nscale)) {
        if (trace >= 1L) message("Computing scaling co-ordinate ", it, ".")
        
        y <- rep(c(-1, 1), length.out=nexamples)
        ynorm.prev <- Inf
        for (k in 1:1000) {
            u <- y / sqrt(y %*% y)
            
            if (nnearest == nexamples) {
                y <- foreach(i=seq_len(nexamples), .combine=c) %dopar% {
                    sum(prox[i, ] * u)
                }
            } else {
                y <- foreach(i=seq_len(nexamples), .combine=c) %dopar% {
                    sum(prox[i, ] * u[prox@examples[i, ]])
                }
            }
            
            y <- 0.5 * (y - (red - sred) * sum(u) - red %*% u)
            
            if (it > 1L) {
                for (j in seq_len(it - 1L)) {
                    y <- y - sum(ev[, j] * u) * dl[j] * ev[, j]
                }
            }
            
            ra <- y %*% u
            sa <- abs(ra)
            ynorm <- sum((y - ra * u) ^ 2)
            
            if (trace >= 2L) message("Scaling co-ordinate ", it, ": ynorm = ",
                                     ynorm, ", threshold = ",
                                     sa * 1.0e-7, ".")
            
            if (ynorm < sa * 1.0e-7 ||
                    (k >= 5L && ynorm / ynorm.prev >= 0.995)) {
                scale[, it] <- sqrt(sa) * u
                ev[, it] <- u
                dl[it] <- ra
				break
            }
            
            ynorm.prev <- ynorm
        }
    }
    
    return(scale)
})
