setGeneric("fastimp", function(forest, ...) standardGeneric("fastimp"))



setMethod("fastimp", signature(forest="bigcforest"), function(forest) {
    # Check argument forest.
    if (!class(forest) == "bigcforest") {
        stop("Argument forest must be a bigcforest created with bigrfc.")
    }
    
    # Calculate fast (gini) importance.
    return(forest@varginidec / mean(forest@varginidec))
})
