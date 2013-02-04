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
    for (t in y@ntrees) {
        x[[x@ntrees + t]] <- y[[1L]]
        y[[1L]] <- NULL
    }
    
    # Set number of trees.
    x@ntrees <- x@ntrees + y@ntrees
    
    # Merge counts of oob times.
    x@oobtimes <- x@oobtimes + y@oobtimes
    
    # Get out-of-bag estimates.
    x@oobvotes <- x@oobvotes + y@oobvotes
    x@avgini <- x@avgini + y@avgini
    
    # Get training set error estimates
    x@oobpred[x@oobtimes > 0L] <-
        sapply(which(x@oobtimes > 0L), function(i) which.max(x@oobvotes[i, ]))
    
    for (c in seq_len(x@ynclass)) {
        x@trainclserr[c] <- sum(class.labels == c & x@oobpred != c)
    }
    
    x@trainerr <- sum(x@trainclserr) / x@nexamples
    x@trainclserr <- x@trainclserr / as.numeric(table(class.labels))
    
    # Calculate confusion matrix
    
    pred <- x@oobpred
    pred[pred == 0L] <- x@ynclass + 1L
    if (length(x@ylevels)) {
        class(pred) <- "factor"
        levels(pred) <- c(x@ylevels, "Never out-of-bag")
        class(class.labels) <- "factor"
        levels(class.labels) <- x@ylevels
    }
    x@trainconfusion <- table(class.labels, pred, dnn=c("Actual", "Predicted"))
    
    return(x)
}

setMethod("merge", signature(x="bigcforest", y="bigcforest"),
          merge.bigcforest)
