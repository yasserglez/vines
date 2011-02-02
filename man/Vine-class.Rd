\name{Vine-class}
\docType{class}

\alias{Vine-class}
\alias{dimnames<-,Vine-method}
\alias{dimnames,Vine-method}
\alias{show,Vine-method}

\title{Base Vine Class}
\description{
  This is the base class of all the classes implementing vine models in
  the package. It is a virtual class, so no objects may be created from it.
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}}
    \item{\code{dimension}:}{Object of class \code{"numeric"}}
    \item{\code{dimensionNames}:}{Object of class \code{"character"}}
    \item{\code{copulas}:}{Object of class \code{"matrix"}}
    \item{\code{trees}:}{Object of class \code{"numeric"}}
  }
}
\section{Methods}{
  \describe{
    \item{dimnames<-}{\code{signature(x = "Vine")}:}
    \item{dimnames}{\code{signature(x = "Vine")}:}
    \item{show}{\code{signature(object = "Vine")}:}
  }
}
\seealso{
  \code{\linkS4class{RVine}}, 
  \code{\linkS4class{CVine}}, 
  \code{\linkS4class{DVine}}.
}
\keyword{classes}
\keyword{multivariate}
\keyword{distribution}
\keyword{models}
