setMethod("merge", signature(x="bigcforest", y="bigcforest"), function(
    x,
    y) {
    
    # Check arguments ----------------------------------------------------------
    
    # Check x and y.
    if (!class(x) == "bigcforest") {
        stop("Argument x must be a bigcforest created with bigrfc.")
    }
    if (!class(y) == "bigcforest") {
        stop("Argument y must be a bigcforest created with bigrfc.")
    }
    # Switch to a more natural naming convention.
    f1 <- x
    f2 <- y
    y <- f1@y
    
    
    
    # Merge forests ------------------------------------------------------------
    
    # Copy trees, keeping memory usage to a minimum.
    oldxntrees <- f1@ntrees
    f1@ntrees <- f1@ntrees + f2@ntrees
    for (t in seq_len(f2@ntrees)) {
        f1[[oldxntrees + t]] <- f2[[1L]]
        f2[[1L]] <- NULL
    }
    
    # Error estimates
    length(f1@trainerr) <- f1@ntrees
    f1@trainclserr <- rbind(f1@trainclserr, matrix(0, f2@ntrees,
                                                   length(levels(f1@y))))
    
    for (t in (oldxntrees + 1):f1@ntrees) {
        tree <- f1[[t]]
        
        # Count number of times examples are out-of-bag.
        w <- tree@insamp == 0L
        f1@oobtimes[w] <- f1@oobtimes[w] + 1L
        
        for (c in seq_along(levels(f1@y))) {
            # Out-of-bag examples with votes for this class.
            w2 <- w & tree@trainpredclass == c
            f1@oobvotes[w2, c] <- f1@oobvotes[w2, c] +
                tree@nodewt[tree@trainprednode[w2]]
        }
        
        f1@oobpred[f1@oobtimes > 0L] <- max.col(f1@oobvotes[f1@oobtimes > 0L, ])
        
        for (c in seq_along(levels(f1@y))) {
            f1@trainclserr[t, c] <- sum(as.integer(y) == c & f1@oobpred != c)
        }
        
        f1@trainerr[t] <- sum(f1@trainclserr[t, ]) / f1@nexamples
        f1@trainclserr[t, ] <- f1@trainclserr[t, ] / as.numeric(f1@ytable)
    }
    
    # Accumulate gini decreases for each variable
    f1@varginidec <- f1@varginidec + f2@varginidec
    
    # Calculate confusion matrix
    
    pred <- f1@oobpred
    pred[pred == 0L] <- length(levels(f1@y)) + 1L
    pred <- factor(pred, levels=seq_len(length(levels(f1@y)) + 1),
                   labels=c(levels(f1@y), "Never out-of-bag"))
    f1@trainconfusion <- table(y, pred, dnn=c("Actual", "Predicted"))
    
    return(f1)
})
