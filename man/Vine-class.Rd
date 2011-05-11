\name{Vine-class}
\docType{class}

\alias{Vine-class}
\alias{show,Vine-method}
\alias{dimnames<-,Vine-method}
\alias{dimnames,Vine-method}

\title{Base Vine Class}

\description{
This is the base class of all the classes implementing vine models in the 
package. It is a virtual class, no objects may be created from it.
}

\section{Slots}{
\describe{
  \item{\code{type}:}{Object of class \code{character}.
    Descriptive name of the type of vine: Canonical vine, D-vine.}
  \item{\code{dimension}:}{Object of class \code{numeric}.
    Dimension of the vine.}
  \item{\code{dimensionNames}:}{Object of class \code{character}.
    Names of the dimensions of the vine. It should be an empty character
    vector or a vector with one element for each dimension of the vine.}
  \item{\code{copulas}:}{Object of class \code{matrix}. This matrix contains
    the \code{\linkS4class{copula}} objects for each copula in the pair-copula
    decomposition. The indexes of this matrix follow the notation used for the 
    subscripts of \eqn{\Theta}{theta} in Aas et al. (2009).}
  \item{\code{trees}:}{Object of class \code{numeric}.
    Number of dependency trees of the vine. It should be an integer between
    \code{0} and \code{dimension - 1} including both endpoints. If \code{trees}
    is less than \code{dimension - 1} the rest of the trees will be assumed 
    to be composed only by Independence copulas.}
}
}

\section{Methods}{
\describe{
  \item{show}{\code{signature(object = "Vine")}: 
    Print a text representation of the vine.}
  \item{dimnames}{\code{signature(x = "Vine")}:
    Retrieve the names of the dimensions of the vine.}
  \item{dimnames<-}{\code{signature(x = "Vine")}:
    Set the names of the dimensions of the vine.}
}
}

\seealso{
\code{\linkS4class{RVine}},
\code{\linkS4class{CVine}},
\code{\linkS4class{DVine}}.
}

\references{
Aas, K. and Czado, C. and Frigessi, A. and Bakken, H. (2009)
Pair-copula constructions of multiple dependence.
\emph{Insurance: Mathematics and Economics} \bold{44}, 182--198.
}

\keyword{classes}
\keyword{multivariate}
\keyword{models}
