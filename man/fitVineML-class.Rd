\name{fitVineML-class}

\docType{class}

\alias{fitVineML-class}
\alias{show,fitVineML-method}

\title{Class for the Results of Vine Inference by Maximum Likelihood}

\description{
Extends the \code{\linkS4class{fitVine}} class to include information
about the Maximum Likelihood inference.
}

\section{Slots}{
See \code{\linkS4class{fitVine}} for information about the inherited slots.

\describe{
  \item{\code{optimMethod}:}{Object of class \code{character}. 
    \code{\link{optim}} method.}
  \item{\code{optimConv}:}{Object of class \code{numeric}. 
    \code{\link{optim}} convergence code. 0 indicates successful completion. 
    If the optimization of the log-likelihood was disabled this will be 0.}
  \item{\code{startParams}:}{Object of class \code{numeric}.
    Vector with the initial parameters of the bivariate copulas in the 
    pair-copula construction calculated using the sequential estimation
    procedure.}
  \item{\code{finalParams}:}{Object of class \code{numeric}.
    Vector with the parameters of the bivariate copulas in the pair-copula 
    construction after the maximization of the log-likelihood. If the 
    optimization of the log-likelihood was disabled this will be equal to 
    the value in the \code{startParams} slot.}
}
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "fitVineML")}:
      Print information about the fitted vine.}
  }
}

\seealso{
  \code{\linkS4class{fitVine}},
  \code{\link{fitVine}},
  \code{\link{vineParameters}}.
}

\keyword{classes}
\keyword{multivariate}
\keyword{models}
