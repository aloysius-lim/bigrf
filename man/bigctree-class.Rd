\name{bigctree-class}
\Rdversion{1.1}
\docType{class}
\alias{bigctree-class}
\title{Classification Trees in Random Forests}
\description{
  Class representing a tree in a classification random forest.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("bigctree", ...)}, but most often are generated when building a \code{\linkS4class{bigcforest}} using \code{\link[bigrf]{bigrfc}}.
}
\section{Slots}{
  \describe{
    \item{\code{insamp}:}{Object of class \code{"integer"}. Number of times each training example is in the bag for this tree.}
    \item{\code{inweight}:}{Object of class \code{"numeric"}. Sum of class weights for each in-bag example.}
    \item{\code{nnodes}:}{Object of class \code{"integer"}. Number of nodes in the tree.}
    \item{\code{treemap}:}{Object of class \code{"matrix"}. Indices of tree nodes. For a given node \code{n}, \code{treemap[n, ]} represents the indices of the left and right children. If a node is terminal, then both elements of \code{treemap[n, ]} are \code{0}.}
    \item{\code{nodeclass}:}{Object of class \code{"integer"}. Majority class of each terminal node.}
    \item{\code{nodewt}:}{Object of class \code{"numeric"}. Sum of class weights of the examples in each terminal node.}
    \item{\code{bestvar}:}{Object of class \code{"integer"}. Best variable that each node was split on, or \code{0} for terminal nodes.}
    \item{\code{bestnumsplit}:}{Object of class \code{"numeric"}. Split point for nodes split on a numerical variable.}
    \item{\code{bestcatsplit}:}{Object of class \code{"list"}. Split point, expressed as integer vectors indicating which factor levels to split on, for nodes split on a categorical variable.}
    \item{\code{termincount}:}{Object of class \code{"numeric"}. Number of training examples in each terminal node.}
    \item{\code{trainprednode}:}{Object of class \code{"integer"}. Terminal node number for each training example.}
    \item{\code{trainpredclass}:}{Object of class \code{"integer"}. Class of terminal node for each training example.}
    \item{\code{tgini}:}{Object of class \code{"numeric"}. Total gini.}
  }
}
\section{Methods}{
  No methods defined with class "bigctree" in the signature.
}
\keyword{classes}
\keyword{models}
\keyword{tree}
