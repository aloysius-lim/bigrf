merge.bigcforest <- function(x, y, class.labels=NULL) {
    # Check arguments ----------------------------------------------------------
    
    # Check x and y.
    if (!class(x) == "bigcforest") {
        stop("Argument x must be a bigcforest created with bigrfc.")
    }
    if (!class(y) == "bigcforest") {
        stop("Argument y must be a bigcforest created with bigrfc.")
    }
    
    # Check class.labels.
    if (x@supervised) {
        if (is.null(class.labels)) {
            stop("Argument class.labels must be specified for supervised ",
                 "learning.")
        }
        if (is.factor(class.labels)) {
            class.labels <- as.integer(class.labels)
        } else if (!is.integer(class.labels)) {
            stop("Argument class.labels must be a factor or integer vector of ",
                 "class labels.")
        }
        if (length(class.labels) != x@nexamples) {
            stop("Argument class.labels must have as many elements as there ",
                 "are rows in x.")
        }
    } else {
        class.labels <- c(rep.int(1L, x@nexamples / 2),
                          rep.int(2L, x@nexamples / 2))
    }
    
    
    
    # Merge forests ------------------------------------------------------------
    
    # Copy trees, keeping memory usage to a minimum.
    oldxntrees <- x@ntrees
    x@ntrees <- x@ntrees + y@ntrees
    for (t in seq_len(y@ntrees)) {
        x[[oldxntrees + t]] <- y[[1L]]
        y[[1L]] <- NULL
    }
    
    # Error estimates
    length(x@trainerr) <- x@ntrees
    x@trainclserr <- rbind(x@trainclserr, matrix(0, y@ntrees, x@ynclass))
    
    for (t in (oldxntrees + 1):x@ntrees) {
        tree <- x[[t]]
        
        # Count number of times examples are out-of-bag.
        w <- tree@insamp == 0L
        x@oobtimes[w] <- x@oobtimes[w] + 1L
        
        for (c in seq_len(x@ynclass)) {
            # Out-of-bag examples with votes for this class.
            w2 <- w & tree@trainpredclass == c
            x@oobvotes[w2, c] <- x@oobvotes[w2, c] +
                tree@nodewt[tree@trainprednode[w2]]
        }
        
        x@oobpred[x@oobtimes > 0L] <- max.col(x@oobvotes[x@oobtimes > 0L, ])
        
        for (c in seq_len(x@ynclass)) {
            x@trainclserr[t, c] <- sum(class.labels == c & x@oobpred != c)
        }
        
        x@trainerr[t] <- sum(x@trainclserr[t, ]) / x@nexamples
        x@trainclserr[t, ] <- x@trainclserr[t, ] / as.numeric(x@ytable)
    }
    
    # Accumulate gini decreases for each variable
    x@varginidec <- x@varginidec + y@varginidec
    
    # Calculate confusion matrix
    
    pred <- x@oobpred
    pred[pred == 0L] <- x@ynclass + 1L
    class(pred) <- "factor"
    levels(pred) <- c(x@ylevels, "Never out-of-bag")
    class(class.labels) <- "factor"
    levels(class.labels) <- x@ylevels
    x@trainconfusion <- table(class.labels, pred, dnn=c("Actual", "Predicted"))
    
    return(x)
}

setMethod("merge", signature(x="bigcforest", y="bigcforest"),
          merge.bigcforest)
