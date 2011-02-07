\name{h-methods}

\docType{methods}

\alias{h}
\alias{h-methods}
\alias{h,copula-method}
\alias{h,indepCopula-method}
\alias{h,normalCopula-method}
\alias{h,tCopula-method}
\alias{h,claytonCopula-method}
\alias{h,gumbelCopula-method}

\title{Methods for the h-functions}

\usage{
h(copula, x, v)
}

\arguments{
\item{copula}{Bivariate \code{\linkS4class{copula}} object.}
\item{x}{\code{numeric} vector.}
\item{v}{\code{numeric} vector.}
}

\description{
  This function should be defined for every copula used in the pair-copula
  constructions and it represents the conditional distribution of \code{x}
  given \code{v} in terms of the copula when \code{x} and \code{v} are Uniform.
  It is defined as the partial derivative of the distribution function of the 
  copula with respect to its second argument.
  \deqn{h(x,v) = F(x|v) = \frac{\partial C(x,v)}{\partial v}}{}
}

\section{Methods}{
\describe{

\item{\code{signature(copula = "copula")}}{
  Default definition of the h-function for a bivariate copula. This method
  is used if there are no definitions specific for a given copula. The
  partial derivative is calculated numerically using the 
  \code{\link{numericDeriv}} function.
}

\item{\code{signature(copula = "indepCopula")}}{
  h-function of the Independence copula.
  \deqn{h(x, v) = x}{}
}

\item{\code{signature(copula = "normalCopula")}}{
  h-function of the Normal copula.
  \deqn{h(x, v; \rho) = 
    \Phi \left( \frac{\Phi^{-1}(x) - \rho\ \Phi^{-1}(v)}
                     {\sqrt{1-\rho^2}} \right)}{}
}

\item{\code{signature(copula = "tCopula")}}{
  h-function of the Student's t copula.
  \deqn{h(x, v; \rho, \nu) = 
    t_{\nu+1} \left( \frac{t^{-1}_{\nu}(x) - \rho\ t^{-1}_{\nu}(v)}
                          {\sqrt{\frac{(\nu+(t^{-1}_{\nu}(v))^2)(1-\rho^2)}
                                      {\nu+1}}} \right)}{}
}

\item{\code{signature(copula = "claytonCopula")}}{
  h-function of the Clayton copula.
  \deqn{h(x, v; \theta) = v^{-\theta-1}(x^{-\theta}+v^{-\theta}-1)^{-1-1/\theta}}{}  
}

\item{\code{signature(copula = "gumbelCopula")}}{
  h-function of the Gumbel copula.
  \deqn{h(x, v; \theta) = C(x, v; \theta)\ \frac{1}{v}\ (-\log v)^{\theta-1}
                          \left((-\log x)^{\theta} + (-\log v)^{\theta} \right)^{1/\theta-1}}{}
}
}}

\keyword{methods}
