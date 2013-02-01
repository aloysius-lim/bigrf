show.bigcforest <- function(object) {
    cat("Random forest with", object@ntrees, "trees, trained on",
        object@nexamples, "examples with", length(object@varselect),
        "variables.\n")
}

setMethod("show", signature(object="bigcforest"), show.bigcforest)
