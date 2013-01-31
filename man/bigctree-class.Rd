\name{bigctree-class}
\Rdversion{1.1}
\docType{class}
\alias{bigctree-class}
\alias{bigctree}
\title{Class \code{"bigctree"}}
\description{
  This class represents a single tree in a random classification forest.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bigctree", ...)}, but most often is the result of a call to \code{\link{bigrfc}}, and stored in a \code{"\linkS4class{bigcforest}"} object.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{insamp}:}{Object of class \code{"integer"} ~~ }
    \item{\code{inweight}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nnodes}:}{Object of class \code{"integer"} ~~ }
    \item{\code{treemap}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{nodeclass}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nodewt}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{bestvar}:}{Object of class \code{"integer"} ~~ }
    \item{\code{bestnumsplit}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{bestcatsplit}:}{Object of class \code{"list"} ~~ }
    \item{\code{termincount}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{trainprednode}:}{Object of class \code{"integer"} ~~ }
    \item{\code{trainpredclass}:}{Object of class \code{"integer"} ~~ }
    \item{\code{tgini}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  No methods defined with class \code{"\linkS4class{bigctree}"} in the signature.
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
  showClass("bigctree")
}
\keyword{classes}
\keyword{models}
\keyword{tree}
