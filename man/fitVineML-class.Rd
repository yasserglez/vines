\name{fitVineML-class}
\docType{class}

\alias{fitVineML-class}
\alias{show,fitVineML-method}

\title{Class for the Results of Vine Inference by Maximum Likelihood}
\description{
  Extends the \code{\linkS4class{fitVine}} class.
}
\section{Slots}{
  See \code{\linkS4class{fitVine}} for information about the inherited slots.

  \describe{
    \item{\code{optimMethod}:}{Object of class \code{"character"} ~~ }
    \item{\code{optimConv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{startParams}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{finalParams}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "fitVineML")}:}
  }
}
\seealso{
  \code{\linkS4class{fitVine}}.
}
\keyword{classes}
\keyword{multivariate}
\keyword{models}
