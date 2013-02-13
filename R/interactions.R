setGeneric("interactions", function(forest, ...) standardGeneric("interactions"))



setMethod("interactions", signature(forest="bigcforest"), function(forest) {
    
    # Check argument forest.
    if (!class(forest) == "bigcforest") {
        stop("Argument forest must be a bigcforest created with bigrfc.")
    }
    
    # Initialize.
    nvar <- length(forest@varselect)
    
    # Rank the variables by gini impurity decrease, for each tree.
    var.rank <- matrix(0, length(forest), nvar)
    for (t in seq_len(length(forest))) {
        tree <- forest[[t]]
        non.zero <- tree@tgini > .Machine$double.eps ^ 0.5
        var.rank[t, non.zero] <- order(tree@tgini[non.zero])
    }
    
    # Histogram of ranks for each variable. Each row represents a variable, and
    # each column represents a rank, where column i is rank i-1 (column 1 is 
    # rank 0, i.e. not ranked).
    hist <- sapply(0:nvar,
                   function(rank)
                       sapply(1:nvar,
                              function(var) sum(var.rank[, var] == rank)))
    hist <- hist / forest@ntrees
    
    # Absolute difference of ranks between pairs of variables.
    effect <- matrix(numeric(), nvar, nvar,
                     dimnames=list(Variable1=names(forest@varselect),
                                   Variable2=names(forest@varselect)))
    for (i in seq_len(nvar)) {
        for (j in i:nvar) {
            effect[i, j] <- sum(abs(var.rank[, i] - var.rank[, j]))
        }
    }
    effect <- effect / forest@ntrees
    
    teffect <- matrix(numeric(), nvar, nvar)
    for (i in seq_len(nvar)) {
        for (j in i:nvar) {
            teffect[i, j] <- 
                sum(sapply(0:nvar,
                           function(rank1)
                               abs(rank1 - 0:nvar) *
                               hist[i, rank1 + 1] * hist[j, ]))
            rcor <- 1 - sum(hist[i, 1:nvar + 1] * hist[j, 1:nvar + 1])
            teffect[i, j] <- teffect[i, j] / rcor
            effect[i, j] <- effect[j, i] <- 100 * (effect[i, j] - teffect[i, j])
        }
    }
    
    return(effect)
})
