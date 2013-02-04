show.bigcforest <- function(object) {
    cat("Random forest with", object@ntrees, "trees, trained on",
        object@nexamples, "examples with", length(object@varselect),
        "variables.\n")
}

setMethod("show", signature(object="bigcforest"), show.bigcforest)



show.bigcprediction <- function(object) {
    cat("Random forest with", object@ntrees, "trees, tested on",
        object@ntest, "examples\n")
    print(object[seq_along(object)])
}

setMethod("show", signature(object="bigcprediction"), show.bigcprediction)
