\name{vineOrdering}

\alias{vineOrdering}

\title{Select an Ordering of the Variables of a Vine}

\description{
Select an ordering of the variables of a vine. The ordering of the
variables determines the bivariate dependences that will be explicit
modeled in the first tree of the vine.
}

\usage{
vineOrdering(type, data, method = "greedy", ...)
}

\arguments{
\item{type}{Type of vine. Supported values are: \code{"CVine"} and 
  \code{"DVine"}.}
\item{data}{Data matrix of pseudo-observations.}
\item{method}{Ordering method. Supported values are: \code{"greedy"}.}
\item{...}{Additional arguments for the ordering method.}
}

\value{
An \code{integer} vector giving the indexes of the variables in order.
This vector should be used to reorder the variables of the \code{data} matrix.
}

\details{
The \code{"greedy"} method uses simple greedy heuristics to select the ordering
of the variables. For the \code{"CVine"} \code{type} it starts by the stronger
pair and then the extreme having the stronger dependence to any other variable 
becomes the first variable in the ordering, which is the variable connected 
to the rest of the variables. For the \code{"DVine"} \code{type} it also starts
by the stronger pair and then iteratively adds variables selecting the variable 
with the stronger dependence to the extremes until the ordering is complete.
Using this heuristics the stronger dependence will be explicit represented in 
the first tree of the vines. Additional parameters for this method are:

\describe{
  \item{\code{according}}{Value used to measure the most dependent pairs.
    The default value is \code{"kendall"} and supported values are:
    \describe{
      \item{\code{"kendall"}}{Absolute value of Kendall's \eqn{\tau}{tau}.}
      \item{\code{"spearman"}}{Absolute value of Spearman's \eqn{\rho}{rho}.}
      \item{\code{"pearson"}}{Absolute value of Pearson's \eqn{\rho}{} correlation coefficient.}
      \item{\code{"df"}}{Smaller degrees of freedom of a bivariate Student's t copula.}
    }
  }
}
}

\examples{
data <- matrix(runif(5 * 100), ncol = 5, nrow = 500)

vineOrdering("CVine", data, method = "greedy", 
             according = "kendall")
vineOrdering("DVine", data, method = "greedy", 
             according = "spearman")
}
