summary.bigcforest <- function(object) {
    cat("Random forest with", object@ntrees, "trees, trained on",
        object@nexamples, "examples with", length(object@varselect),
        "variables.\n\n")
    cat("Training set labels:\n")
    print(object@ytable)
    
    cat("\nOverall error rate:",
        format(100 * object@trainerr, digits=3, nsmall=2), "\n\n")
    
    cat("Training set confusion matrix (OOB):\n")
    print(object@trainconfusion)
}

setMethod("summary", signature(object="bigcforest"), summary.bigcforest)
