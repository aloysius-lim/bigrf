show.bigcprediction <- function(object) {
    cat("Random forest with", object@ntrees, "trees, tested on",
        object@ntest, "samples.\n")
    print(object[seq_along(object)])
}

setMethod("show", signature(object="bigcprediction"), show.bigcprediction)
