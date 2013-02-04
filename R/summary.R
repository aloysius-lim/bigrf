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



summary.bigcprediction <- function(object) {
    cat("Random forest with", object@ntrees, "trees, tested on",
        object@ntest, "examples\n\n")
    if (object@testlabelled) {
        cat("Test set labels:\n")
        print(object@testytable)
        
        cat("\nOverall error rate:",
            format(100 * object@testerr, digits=3, nsmall=2), "\n\n")
        
        cat("Test set confusion matrix (OOB):\n")
        print(object@testconfusion)
    }
}

setMethod("summary", signature(object="bigcprediction"), summary.bigcprediction)
