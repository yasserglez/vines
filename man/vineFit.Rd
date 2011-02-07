\name{vineFit}

\alias{vineFit}

\title{Vine Inference}

\description{
Fit a vine model to multivariate data in to the unit hypercube. The 
data can be pseudo-observations constructed from empirical or parametric
marginal cumulative distribution functions.
}

\usage{
vineFit(type, data, method = "ml", ...)
}

\arguments{
\item{type}{Type of vine. Supported values: are \code{"CVine"} 
  and \code{"DVine"}.}
\item{data}{Data matrix of pseudo-observations.}
\item{method}{Inference method. Supported values are: 
  \code{"ml"} (Maximum Likelihood).}
\item{...}{Additional arguments for the inference method.}
}

\value{
An object of class \code{\linkS4class{vineFit}} or a subclass with 
information specific to the used inference method. The \code{vine} slot 
of this object will contain the fitted \code{\linkS4class{Vine}} object.
}

\details{

The \code{"ml"} (Maximum Likelihood) method uses sequential estimation 
procedure outlined in Aas et al. (2009) to determine the parametric
shape and initial values for the parameters of each bivariate copula 
and then execute a numerical maximization of the full log-likelihood 
of the vine. Additional arguments for this method are:

\describe{
  \item{\code{selectCopula}}{Function called to select the bivariate copulas
    in the pair-copula construction. This function receives a bivariate
    sample and should return a \code{\linkS4class{copula}} object that fits
    this sample. The default value is \code{function (x, y) indepCopula()}
    meaning to use Independence copulas for all the bivariate copulas.}
  \item{\code{trees}}{Maximum number of dependency trees of the vine. The
    \code{selectCopula} function will be executed for copulas in the first
    \code{trees} trees of the vine and Independence copulas will be used
    for the rest of the trees. The default value is \code{ncol(data) - 1}
    meaning to capture dependence in all the trees of the vine.}
  \item{\code{truncMethod}}{Method to automatically truncate the vine if enough
    dependence is capture in a given number of trees. Supported methods are 
    \code{"AIC"} and \code{"BIC"}. See Brechmann et al. (2010) for information 
    about these methods. The default value is \code{""} meaning no truncation.}
  \item{\code{optimMethod}}{\code{\link{optim}} method used in the optimization
    of the log-likelihood function. If \code{""} is specified the optimization
    is disabled and the vine with the initial parameters calculated using the
    sequential estimation procedure is returned. The default value is
    \code{"Nelder-Mead"}.}
  \item{\code{optimControl}}{List of control parameters for \code{\link{optim}}.
    The default value is \code{list()}.}
}
}

\examples{
data <- matrix(runif(5 * 100), ncol = 5, nrow = 100)
colnames(data) <- c("A", "B", "C", "D", "E")

selectCopula <- function (x, y) {
    data <- cbind(x, y)
    fit <- fitCopula(normalCopula(0), data, method = "itau")
    fit@copula
}
fit <- vineFit("DVine", data, method = "ml", 
               selectCopula = selectCopula, 
               optimMethod = "")

show(fit)
show(fit@vine)
}

\seealso{
\code{\linkS4class{CVine}},
\code{\linkS4class{DVine}},
\code{\linkS4class{vineFit}},
\code{\linkS4class{vineFitML}}.
}

\references{
Aas, K. and Czado, C. and Frigessi, A. and Bakken, H. (2009)
Pair-copula constructions of multiple dependence.
\emph{Insurance: Mathematics and Economics} \bold{44}, 182--198.

Brechmann, E. C. and Czado, C. and Aas, K. (2010) Truncated regular 
vines in high dimensions with application to financial data. 
\emph{Norwegian Computing Center, NR}. Note SAMBA/60/10.
}
