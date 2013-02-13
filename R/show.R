setMethod("show", signature(object="bigcforest"), function(object) {
    cat("Random forest with", object@ntrees, "trees, trained on",
        object@nexamples, "examples with", length(object@varselect),
        "variables.\n")
    
    cat("Cumulative error rates:\n")
    errors <- cbind(object@trainerr, object@trainclserr)
    dimnames(errors)[[2]][1] <- "Overall"
    print(errors)
})



setMethod("show", signature(object="bigcprediction"), function(object) {
    cat("Predictions on", object@ntest, "examples using random forest with",
        object@ntrees, "trees.\n")
    print(object[seq_along(object)])
})
