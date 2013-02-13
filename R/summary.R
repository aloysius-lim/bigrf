setMethod("summary", signature(object="bigcforest"), function(object) {
    cat("Random forest with", object@ntrees, "trees, trained on",
        object@nexamples, "examples with", length(object@varselect),
        "variables.\n\n")
    
    cat("Training set labels:\n")
    print(object@ytable)
    cat("\n")
    
    cat("Training class weights:\n")
    print(object@yclasswts)
    cat("\n")
    
    cat("Error rates:\n")
    errors <- c(object@trainerr[object@ntrees],
                object@trainclserr[object@ntrees, ]) * 100
    names(errors) <- c("Overall", levels(object@y))
    print(errors)
    cat("\n")
    
    cat("Training set confusion matrix (OOB):\n")
    print(object@trainconfusion)
})



setMethod("summary", signature(object="bigcprediction"), function(object) {
    cat("Predictions on", object@ntest, "examples using random forest with",
        object@ntrees, "trees.\n\n")
    
    if (object@testlabelled) {
        cat("Test set labels:\n")
        print(object@testytable)
        
        cat("\nOverall error rate:",
            format(100 * object@testerr, digits=3, nsmall=2), "\n\n")
        
        cat("Test set confusion matrix (OOB):\n")
        print(object@testconfusion)
    } else {
        cat("Test set was not labelled.\n")
    }
})