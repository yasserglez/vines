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

setGeneric("hinverse", 
    function (copula, u, v) {
      if (copula@dimension != 2) {
        stop("inverse of h-functions only defined for bivariate copulas")
      }
      # Avoid numerical problems at the boundary of the interval.
      eps <- .Machine$double.eps^0.5
      u[u < eps] <- eps
      u[(1 - u) < eps] <- 1 - eps
      v[v < eps] <- eps
      v[(1 - v) < eps] <- 1 - eps
      standardGeneric("hinverse")
    },
    signature = "copula")


# For the expression for the Gaussian, Student's t, Clayton and Gumbel copulas
# see Aas, K., Czado, C., Frigessi, A. & Bakken, H. Pair-copula constructions 
# of multiple dependence. Insurance Mathematics and Economics, 2007, Vol. 44, 
# pp. 182-198.


hinverseCopula <- function (copula, u, v) {
  # Last resort is to evaluate the inverse of the h-function numerically.
  f <- function (x, u, v) h(copula, x, v) - u
  root <- function (u, v) {
    eps <- .Machine$double.eps^0.5
    if (u == eps || u == 1 - eps) u
    else {
      # Choose one of the roots in the interval.
      sample(uniroot.all(f, lower = eps, upper = 1 - eps, u = u, v = v), 1)
    }
  }
  sapply(seq(along = u), function (i) root(u[i], v[i]))
}

setMethod("hinverse", "copula", hinverseCopula)


hinverseIndepCopula <- function (copula, u, v) {
  u
}

setMethod("hinverse", "indepCopula", hinverseIndepCopula)


hinverseNormalCopula <- function (copula, u, v) {
  rho <- copula@parameters
  pnorm(qnorm(u) * sqrt(1 - rho^2) + rho*qnorm(v))
}

setMethod("hinverse", "normalCopula", hinverseNormalCopula)


hinverseClaytonCopula <- function (copula, u, v) {
  theta <- copula@parameters
  ((u * v^(theta+1)) ^ (-theta/(theta+1)) + 1 - v^(-theta)) ^ (-1/theta)
}

setMethod("hinverse", "claytonCopula", hinverseClaytonCopula)


hinverseTCopula <- function (copula, u, v) {
  rho <- copula@parameters[1]
  df <- if (copula@df.fixed) copula@df else copula@parameters[2]
  q <- qt(u, df+1) * 
      sqrt(((df + qt(v, df)^2) * (1 - rho^2)) / (df+1)) + 
      rho*qt(v, df)
  pt(q, df)
}

setMethod("hinverse", "tCopula", hinverseTCopula)
