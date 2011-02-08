\name{hinverse-methods}

\docType{methods}

\alias{hinverse}
\alias{hinverse-methods}
\alias{hinverse,copula-method}
\alias{hinverse,indepCopula-method}
\alias{hinverse,normalCopula-method}
\alias{hinverse,tCopula-method}
\alias{hinverse,claytonCopula-method}

\title{Methods for the Inverse of the h-functions}

\usage{
hinverse(copula, u, v)
}

\arguments{
\item{copula}{A bivariate \code{\linkS4class{copula}} object.}
\item{u}{\code{numeric} vector.}
\item{v}{\code{numeric} vector.}
}

\description{
  This function should be defined for every copula used in the pair-copula
  constructions or evaluated numerically. It represents the inverse of the 
  h-function, or equivalently the inverse of the conditional distribution 
  function.
}

\section{Methods}{
\describe{

\item{\code{signature(copula = "copula")}}{
  Default definition of the inverse of the h-function for a bivariate copula. 
  This method is used if there are no definitions specific for a given copula. 
  The inverse is calculated numerically using the \code{\link{uniroot}}
  function.
}

\item{\code{signature(copula = "indepCopula")}}{
  Inverse of the h-function of the Independence copula.
  \deqn{h^{-1}(u, v) = u}{}
}

\item{\code{signature(copula = "normalCopula")}}{
  Inverse of the h-function of the Normal copula.
  \deqn{h^{-1}(u, v; \rho) = \Phi \left( \Phi^{-1}(u) \sqrt{1-\rho^2} +
                                         \rho\ \Phi^{-1}(v) \right)}{}  
}

\item{\code{signature(copula = "tCopula")}}{
  Inverse of the h-function of the Student's t copula.
  \deqn{h^{-1}(u, v; \rho, \nu) = 
    t_{\nu} \left( t^{-1}_{\nu+1}(u)\
                   \sqrt{\frac{(\nu+(t^{-1}_{\nu}(v))^2)(1-\rho^2)}{\nu+1}} +
                   \rho\ t^{-1}_{\nu}(v)
            \right)}{}
}

\item{\code{signature(copula = "claytonCopula")}}{
  Inverse of the h-function of the Clayton copula.
  \deqn{h^{-1}(u, v; \theta) = \left( \left( u\ v^{\theta+1}\right)^{-\frac{\theta}{\theta+1}} + 
                                      1 - v^{-\theta} \right)^{-1/\theta}}{}  
}
}}

\keyword{methods}
