setGeneric("proximities", function(forest, ...) standardGeneric("proximities"))



setMethod("proximities", signature(forest="bigcforest"), function(
    forest,
    nprox) {
    
    prox <- big.matrix(forest@nexamples, nprox)
    
})
