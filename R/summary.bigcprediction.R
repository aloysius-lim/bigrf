summary.bigcprediction <- function(object) {
    cat("Random forest with", object@ntrees, "trees, tested on",
        object@ntest, "samples.\n\n")
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
