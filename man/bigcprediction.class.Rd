\name{bigcprediction-class}
\Rdversion{1.1}
\docType{class}
\alias{bigcprediction-class}
\alias{show,bigcprediction-method}
\alias{summary,bigcprediction-method}

\title{Class \code{"bigcprediction"}}
\description{
The outputs of \code{\link[bigrf]{predict}} on a \code{"\linkS4class{bigcforest}"}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bigcprediction", ...)}, but most often is the result of a call to \code{\link[bigrf]{predict}}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"integer"}. The predicted class for each test example.}
    \item{\code{ntest}:}{Object of class \code{"integer"}. The number of test examples used for prediction.}
    \item{\code{testlabelled}:}{Object of class \code{"logical"}. Whether the test examples were labelled. If TRUE, then error estimates and the confusion matrix are available.}
    \item{\code{ynclass}:}{Object of class \code{"integer"}. Number of classes.}
    \item{\code{ntrees}:}{Object of class \code{"integer"}. Number of trees in the random forest.}
    \item{\code{testytable}:}{Object of class \code{"table"}. Counts of test examples in each class, if test examples were labelled. Otherwise, \code{NULL}.}
    \item{\code{testvotes}:}{Object of class \code{"matrix"}. Weighted class votes for each test example. The prediction for each example is the class that received the highest total vote.}
    \item{\code{testclserr}:}{Object of class \code{"numeric"}. Prediction error for each class, if test examples were labelled. Otherwise, \code{NULL}.}
    \item{\code{testerr}:}{Object of class \code{"numeric"}. Total prediction error for all classes, if test examples were labelled. Otherwise, \code{NULL}.}
    \item{\code{testconfusion}:}{Object of class \code{"table"}. The confusion matrix for the test set, if test examples were labelled. Otherwise, \code{NULL}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{integer}"}, from data part.
Class \code{"\linkS4class{numeric}"}, by class "integer", distance 2.
Class \code{"\linkS4class{vector}"}, by class "integer", distance 2.
Class \code{"\linkS4class{data.frameRowLabels}"}, by class "integer", distance 2.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "bigcprediction")}: Prints the prediction results.}
    \item{summary}{\code{signature(object = "bigcprediction")}: Prints summary information on the prediction results, including test error estimates and the confusion matrix, if test labells were supplied during prediction.}
	 }
}
\examples{
showClass("bigcprediction")
}
\keyword{classes}
