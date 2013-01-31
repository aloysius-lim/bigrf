\name{show-methods}
\docType{methods}
\alias{show,bigcforest-method}
\title{Methods for Function \code{show} in Package \pkg{bigrf}}
\description{
  Prints information about a random forest.
}
\usage{
  \S4method{show}{bigcforest}(object)
}
\arguments{
  \item{object}{
    A random forest of class \code{"\linkS4class{bigcforest}"}.
  }
}
\section{Methods}{
\describe{
\item{\code{signature(object = "bigcforest")}}{
Prints information about the given random forest, displaying summary error statistics and confusion matrices.
}
}}
\keyword{methods}
\keyword{print}
