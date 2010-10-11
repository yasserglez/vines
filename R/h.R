# vines: GNU R package for multivariate dependence modeling with vines
# Copyright (C) 2010 Yasser González Fernández <ygonzalezfernandez@gmail.com>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details.
#
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>.

setGeneric("h", 
    function (copula, x, v) {
      if (copula@dimension != 2) {
        stop("h-functions only defined for bivariate copulas")
      }
      # Avoid numerical problems at the boundary of the interval.
      eps <- .Machine$double.neg.eps
      x[x < eps] <- eps
      x[x > 1 - eps] <- 1 - eps
      v[v < eps] <- eps
      v[v > 1 - eps] <- 1 - eps
      r <- standardGeneric("h")
      r[r < eps] <- eps
      r[r > 1 - eps] <- 1 - eps
      r
    },
    signature = "copula")


# For the expression for the Gaussian, Student's t, Clayton and Gumbel copulas
# see Aas, K., Czado, C., Frigessi, A. & Bakken, H. Pair-copula constructions 
# of multiple dependence. Insurance Mathematics and Economics, 2007, Vol. 44, 
# pp. 182-198.


hCopula <- function (copula, x, v) {
  # Last resort is to evaluate the h-function numerically.
  e <- new.env()
  assign("x", x, envir = e)
  assign("v", v, envir = e)
  r <- numericDeriv(quote(pcopula(copula, cbind(x, v))), c("v"), e)
  diag(attr(r, "gradient"))
}

setMethod("h", "copula", hCopula)


hIndepCopula <- function (copula, x, v) {
  x
}

setMethod("h", "indepCopula", hIndepCopula)


hNormalCopula <- function (copula, x, v) {
  rho <- copula@parameters
  pnorm((qnorm(x) - rho*qnorm(v)) / sqrt(1 - rho^2))
}

setMethod("h", "normalCopula", hNormalCopula)


hClaytonCopula <- function (copula, x, v) {
  theta <- copula@parameters
  v^(-theta-1) * (x^(-theta) + v^(-theta) - 1) ^ (-1 - 1/theta)
}

setMethod("h", "claytonCopula", hClaytonCopula)


hGumbelCopula <- function (copula, x, v) {
  theta <- copula@parameters
  pcopula(copula, cbind(x, v)) * 1/v * (-log(v)) ^ (theta-1) * 
      ((-log(x))^theta + (-log(v))^theta) ^ (1/theta-1)
}

setMethod("h", "gumbelCopula", hGumbelCopula)


hTCopula <- function (copula, x, v) {
  rho <- copula@parameters[1]
  df <- if (copula@df.fixed) copula@df else copula@parameters[2]
  q <- (qt(x, df) - rho*qt(v, df)) / 
      sqrt(((df + qt(v, df)^2) * (1 - rho^2)) / (df+1))
  pt(q, df+1)
}

setMethod("h", "tCopula", hTCopula)
