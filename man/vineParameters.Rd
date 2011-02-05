\name{vineParameters}

\alias{vineParameters}
\alias{vineParameters<-}

\title{Parameters of a Vine}

\description{
Retrieve or set the parameters of a \code{\linkS4class{Vine}}.
}

\usage{
parameters <- vineParameters(vine)
vineParameters(vine) <- parameters
}

\arguments{
\item{vine}{A \code{\linkS4class{Vine}} object.}
\item{parameters}{\code{numeric} vector with the parameters of all the 
  copulas in the pair-copula construction. This vector is the result
  of the concatenation of the \code{parameters} slot of the 
  \code{\linkS4class{copula}} objects in the \code{copulas} slot of the 
  \code{\linkS4class{Vine}} object (by columns).}
}

\examples{
dimension <- 3
copulas <- matrix(list(gumbelCopula(3.45), claytonCopula(2.75),
                       tCopula(0.75, df = 5), NULL),
                  ncol = dimension - 1, nrow = dimension - 1)
vine <- DVine(dimension = dimension, trees = dimension - 1,
              copulas = copulas)
show(vine)
vineParameters(vine) <- vineParameters(vine) / 2
show(vine) 
}
