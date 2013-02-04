setGeneric("fastimp", function(forest, ...) standardGeneric("fastimp"))



fastimp.bigcforest <- function(forest) {
    # Check arguments ----------------------------------------------------------
    
    # Check forest.
    if (!class(forest) == "bigcforest") {
        stop("Argument forest must be a bigcforest created with bigrfc.")
    }
    
    
    
    # Calculate fast (gini) importance -----------------------------------------
    return(forest@varginidec / mean(forest@varginidec))
}

setMethod("fastimp", signature(forest="bigcforest"),  fastimp.bigcforest)
