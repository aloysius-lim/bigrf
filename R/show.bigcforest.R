show.bigcforest <- function(object) {
    cat("Random forest with", object@ntrees, "trees, trained on",
        object@nsample, "samples with", length(object@varselect),
        "variables.\n\n")
    cat("Overall error rate:",
        format(100 * object@trainerr, digits=3, nsmall=2), "\n\n")
    
    # if (ntest > 0 && labelts != 0) {
    #   write[*,*] 'final error test %    ',100*errts 
    # } 
    
    cat("Training set confusion matrix (OOB):\n")
    print(object@trainconfusion)
    
    # if (ntest > 0 && labelts != 0) {
    #   mtab <- table(clts, jests, dnn=c("Actual", "Predicted"))
    #   cat("Test set confusion matrix:\n")
    #   print(mtab)
    # }   
}

setMethod("show", signature(object="bigcforest"), show.bigcforest)
