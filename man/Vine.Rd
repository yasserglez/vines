\name{Vine}

\alias{Vine}
\alias{dvine}
\alias{dvine-methods}
\alias{dvine,CVine-method}
\alias{dvine,DVine-method}
\alias{pvine}
\alias{pvine-methods}
\alias{pvine,CVine-method}
\alias{pvine,DVine-method}
\alias{rvine}
\alias{rvine-methods}
\alias{rvine,CVine-method}
\alias{rvine,DVine-method}

\title{Vine Distribution Functions}

\description{
   Density evaluation, distribution function evaluation, and random generation 
   for vines.
}

\usage{
  dvine(vine, u)
  pvine(vine, u)
  rvine(vine, n)
}

\arguments{
  \item{vine}{A \code{\linkS4class{Vine}} object.}
  \item{u}{A vector of the same dimension of the vine or a matrix with one 
    column for each dimension of the vine and rows giving the coordinates 
    of the points where the density or distribution functions needs to be 
    evaluated.}
  \item{n}{Number of observations to be generated.}
}

\details{
  All this functions are S4 generic functions on the first argument.

  The implementation of the \code{dvine} function for \code{\linkS4class{CVine}} 
  and \code{\linkS4class{DVine}} is based on the Algorithms 3 and 4 from
  Aas et al. (2009).
  
  The \code{pvine} function is currently evaluated through the numerical 
  integration of the density function. This is slow even for small dimensions.
  See the documentation of the \code{\link{cubature}} package for information
  about the algorithm used.
  
  The implementation of the \code{rvine} function for \code{\linkS4class{CVine}} 
  and \code{\linkS4class{DVine}} is based on the Algorithms 1 and 2 from
  Aas et al. (2009).
}

\value{
  \code{dnorm} returns a vector with the evaluation of the density. 
  \code{pnorm} returns a vector with the evaluation of the distribution 
  function. \code{rvine} returns a matrix with one column for each dimension
  of the vine and each row corresponding to a simulation of the vine.
}

\references{
  Aas, K. and Czado, C. and Frigessi, A. and Bakken, H. (2009)
  Pair-copula constructions of multiple dependence.
  \emph{Insurance: Mathematics and Economics} \bold{44}, 182--198.

  Bedford, T. and Cooke, R. M. (2001) Monte Carlo simulation of vine
  dependent random variables for applications in uncertainty analysis. 
  In \emph{2001 Proceedings of ESREL2001, Turin, Italy}.

  Bedford, T. and Cooke, R. M. (2001) Probability density decomposition 
  for conditionally dependent random variables modeled by vines. 
  \emph{Annals of Mathematics and Artificial Intelligence} \bold{32},
  245--268.

  Kurowicka, D. and Cooke, R. M. (2005) Sampling algorithms for generating 
  joint uniform distributions using the vine-copula method. In \emph{3rd 
  IASC world conference on Computational Statistics & Data Analysis, 
  Limassol, Cyprus}.
}

\examples{
copulas <- matrix(list(normalCopula(0.5), claytonCopula(2.75),
                       tCopula(0.75, df = 2), NULL), 2, 2)
vine <- DVine(dimension = 3, trees = 2, copulas = copulas)
dimnames(vine) <- c("A", "B", "C")
show(vine)
data <- rvine(vine, 1)
dvine(vine, data)
pvine(vine, data)
}

\keyword{methods}
\keyword{distribution}
\keyword{multivariate}