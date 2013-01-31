\name{bigcforest-class}
\Rdversion{1.1}
\docType{class}
\alias{bigcforest-class}
\alias{bigcforest}
\title{Class \code{"bigcforest"}}
\description{
  This class represents a random classification forest.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bigcforest", ...)}, but most often is the result of a call to \code{\link{bigrfc}}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
    \item{\code{supervised}:}{Object of class \code{"logical"} ~~ }
    \item{\code{nsample}:}{Object of class \code{"integer"} ~~ }
    \item{\code{factors}:}{Object of class \code{"logical"} ~~ }
    \item{\code{ylevels}:}{Object of class \code{"character"} ~~ }
    \item{\code{nlevels}:}{Object of class \code{"integer"} ~~ }
    \item{\code{ytable}:}{Object of class \code{"table"} ~~ }
    \item{\code{varselect}:}{Object of class \code{"integer"} ~~ }
    \item{\code{contvarseq}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nclass}:}{Object of class \code{"integer"} ~~ }
    \item{\code{classweights}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{ntrees}:}{Object of class \code{"integer"} ~~ }
    \item{\code{maxnodes}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nsplitvar}:}{Object of class \code{"integer"} ~~ }
    \item{\code{maxndsize}:}{Object of class \code{"integer"} ~~ }
    \item{\code{maxeslevels}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nrandsplit}:}{Object of class \code{"integer"} ~~ }
    \item{\code{oobtimes}:}{Object of class \code{"integer"} ~~ }
    \item{\code{oobvotes}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{oobpred}:}{Object of class \code{"integer"} ~~ }
    \item{\code{trainclserr}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{trainerr}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{trainconfusion}:}{Object of class \code{"table"} ~~ }
    \item{\code{avgini}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{printerrfreq}:}{Object of class \code{"integer"} ~~ }
    \item{\code{printclserr}:}{Object of class \code{"logical"} ~~ }
    \item{\code{cachepath}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{list}"}, from data part.
  Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{grow}{\code{signature(forest = "bigcforest")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
  showClass("bigcforest")
}
\keyword{classes}
\keyword{models}
\keyword{tree}

