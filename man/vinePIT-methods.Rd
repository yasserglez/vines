\name{vinePIT-methods}
\docType{methods}

\alias{vinePIT-methods}
\alias{vinePIT,CVine-method}
\alias{vinePIT,DVine-method}

\title{Vine Probability Integral Transform Methods}

\description{
Probability Integral Transform (PIT) of Rosenblatt (1952) for vine models. 
The PIT converts a set of dependent variables into a new set of variables 
which are independent and Uniform in \eqn{(0,1)}{(0,1)} under the hypothesis 
that the data originate from a given multivariate distribution.
}

\usage{
vinePIT(vine, u)
}

\arguments{
\item{vine}{A \code{\linkS4class{Vine}} object.}
\item{u}{A vector of the same dimension of the vine or a matrix with one 
  column for each dimension of the vine.}
}

\value{
A matrix with one column for each dimension of the vine and each row 
corresponding to a transformed observation.
}

\section{Methods}{
\describe{
  \item{\code{signature(vine = "CVine")}}{
    Implementation of the PIT algorithm for \code{\linkS4class{CVine}} objects. 
    Based on the Algorithm 5 from Aas et al. (2009).
  }

  \item{\code{signature(vine = "DVine")}}{
    Implementation of the PIT algorithm for \code{\linkS4class{DVine}} objects. 
    Based on the Algorithm 6 from Aas et al. (2009).
  }
}}

\references{
Aas, K. and Czado, C. and Frigessi, A. and Bakken, H. (2009)
Pair-copula constructions of multiple dependence.
\emph{Insurance: Mathematics and Economics} \bold{44}, 182--198.

Rosenblatt, M. (1952) Remarks on multivariate transformation.
\emph{Annals of Mathematical Statistics} \bold{23}, 1052--1057.
}

\examples{
dimension <- 3
copulas <- matrix(list(normalCopula(0.5), 
                       claytonCopula(2.75),
                       NULL, NULL),
                  ncol = dimension - 1,
                  nrow = dimension - 1)
vine <- CVine(dimension = dimension, trees = 1,
              copulas = copulas)

data <- matrix(runif(dimension * 100), 
               ncol = dimension, nrow = 100)

vinePIT(vine, data)
}

\keyword{methods}
