\name{vineOrder}

\alias{vineOrder}

\title{Select an Order of the Variables of a Vine}

\description{
Select an order of the variables of a vine. The order of the
variables determines the bivariate dependences that will be explicit
modeled in the first tree of the vine.
}

\usage{
vineOrder(type, data, method = "greedy", ...)
}

\arguments{
\item{type}{Type of vine. Supported values are: 
  \code{"CVine"} and \code{"DVine"}.}
\item{data}{Data matrix of pseudo-observations.}
\item{method}{Ordering method. Supported values are: 
  \code{"random"}, \code{"greedy"}.}
\item{...}{Additional arguments for the ordering method.}
}

\value{
An \code{integer} vector with the indexes of the variables in order.
This vector can be used to reorder the variables of the \code{data} matrix.
}

\details{
The order of the variables of the first tree of a \code{"DVine"} determines 
the structure of all the trees of the vine. This is not the case for the 
\code{"CVine"} where you can select the root node of each tree in the vine, 
but this is currently not implemented in the package.

The \code{"random"} method returns a random permutation of the variables.

The \code{"greedy"} method returns an order of the variables that tries to
capture as much dependency as possible in the first tree of the vine. This
is done by finding the order of the variables that defines the tree that
maximizes a given dependency measure used as edge weights. For the 
\code{"CVine"} this is done iteratively checking each variable as the 
root node. For the \code{"DVine"} this is equivalent to solve the 
\emph{Traveling Salesman Problem} (TSP) and it is done using the 
\emph{Cheapest Insertion Algorithm} implemented by the \code{\link{solve_TSP}} 
function from the \pkg{TSP} package. See Brechmann, E. C. (2010) for 
details. Additional parameters for this method are:

\describe{
  \item{\code{according}}{Dependence measure. The default value is 
    \code{"kendall"} and supported values are:
    \describe{
      \item{\code{"kendall"}}{Absolute value of Kendall's \eqn{\tau}{tau}.}
      \item{\code{"spearman"}}{Absolute value of Spearman's \eqn{\rho}{rho}.}
      \item{\code{"pearson"}}{Absolute value of Pearson's product moment
        correlation coefficient.}
      \item{\code{"df"}}{Smaller degrees of freedom of a bivariate
        Student's t copula.}
    }
  }
}
}

\examples{
data <- matrix(runif(5 * 100), ncol = 5, nrow = 500)

vineOrder("CVine", data, method = "random")
vineOrder("DVine", data, method = "greedy",
          according = "spearman")
}

\references{
Brechmann, E. C. (2010) Truncated and simplified regular vines 
and their applications. Diploma Thesis. \emph{Technische 
Universität München}.
}
