setMethod("show", signature(object="bigcforest"), function(object) {
    cat("Random forest with", object@ntrees, "trees, trained on",
        object@nexamples, "examples with", length(object@varselect),
        "variables.\n")
})



setMethod("show", signature(object="bigcprediction"), function(object) {
    cat("Random forest with", object@ntrees, "trees, tested on",
        object@ntest, "examples\n")
    print(object[seq_along(object)])
})
