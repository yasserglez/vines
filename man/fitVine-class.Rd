\name{fitVine-class}

\docType{class}

\alias{fitVine-class}
\alias{show,fitVine-method}

\title{Class for the Results of Vine Inference}

\description{
This is the base class of all classes giving information about vine inference.
Objects of this class (or subclasses containing information about specific
inference methods) are created as results of calls to the \code{\link{fitVine}}
function.
}

\section{Slots}{
\describe{
  \item{\code{vine}:}{Object of class \code{\linkS4class{Vine}}.
    Fitted vine.}
  \item{\code{observations}:}{Object of class \code{numeric}.
    Number of observations in the sample.}
  \item{\code{method}:}{Object of class \code{character}.
    Inference method.}
}
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "fitVine")}: 
      Print information about the fitted vine.}
  }
}

\seealso{
  \code{\link{fitVine}},
  \code{\linkS4class{fitVineML}}.
}

\keyword{classes}
\keyword{multivariate}
\keyword{models}
