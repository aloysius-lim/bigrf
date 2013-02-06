setMethod("summary", signature(object="bigcforest"), function(object) {
    cat("Random forest with", object@ntrees, "trees, trained on",
        object@nexamples, "examples with", length(object@varselect),
        "variables.\n\n")
    cat("Training set labels:\n")
    print(object@ytable)
    cat("\n")
    
    cat("Error rates:\n")
    errors <- c(object@trainerr[object@ntrees],
                object@trainclserr[object@ntrees, ]) * 100
    names(errors) <- c("Overall", object@ylevels)
    print(errors)
    cat("\n")
    
    cat("Training set confusion matrix (OOB):\n")
    print(object@trainconfusion)
})


setMethod("summary", signature(object="bigcprediction"), function(object) {
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
})