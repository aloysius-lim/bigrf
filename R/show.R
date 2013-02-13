setMethod("show", signature(object="bigcforest"), function(object) {
    cat("Random forest with", object@ntrees, "trees, trained on",
        object@nexamples, "examples with", length(object@varselect),
        "variables.\n")
})



setMethod("show", signature(object="bigcprediction"), function(object) {
    cat("Predictions on ", object@ntest, " examples using random forest with",
        object@ntrees, ".\n")
    print(object[seq_along(object)])
})
