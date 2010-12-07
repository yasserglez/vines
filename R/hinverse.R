# vines: R package for multivariate dependence modeling with vines
# Copyright (C) 2010 Yasser González-Fernández <ygonzalezfernandez@gmail.com>
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
    function (copula, u, v) standardGeneric("hinverse"),
    signature = "copula")


hinverseCopula <- function (copula, u, v) {
  # Last resort is to evaluate the inverse of the h-function numerically.
  zero <- .Machine$double.eps
  one <- 1 - .Machine$double.neg.eps

  u[u < zero] <- zero
  u[u > one] <- one
  v[v < zero] <- zero
  v[v > one] <- one

  f <- function (x, u, v) abs(h(copula, x, v) - u)
  z <- function (i) optimize(f, c(zero, one), tol = 0.01, u = u[i], v = v[i])$minimum
  r <- sapply(seq(along = u), z)

  r
}

setMethod("hinverse", "copula", hinverseCopula)


hinverseIndepCopula <- function (copula, u, v) {
  u
}

setMethod("hinverse", "indepCopula", hinverseIndepCopula)


# See Aas, K., Czado, C., Frigessi, A. & Bakken, H. Pair-copula constructions 
# of multiple dependence. Insurance Mathematics and Economics, 2009, Vol. 44, 
# pp. 182-198 for the expression for the Gaussian, Student's t, Clayton and 
# Gumbel copulas.

hinverseNormalCopula <- function (copula, u, v) {
  zero <- .Machine$double.eps
  one <- 1 - .Machine$double.neg.eps
  
  u[u < zero] <- zero
  u[u > one] <- one
  v[v < zero] <- zero
  v[v > one] <- one
  
  rho <- copula@parameters
  rho[rho == -1] <- -1 + .Machine$double.eps
  rho[rho == 1] <- 1 - .Machine$double.neg.eps
  
  r <- pnorm(qnorm(u) * sqrt(1 - rho^2) + rho*qnorm(v))
  
  r[u <= zero | r < zero] <- zero
  r[u >= one | r > one] <- one

  r
}

setMethod("hinverse", "normalCopula", hinverseNormalCopula)


hinverseTCopula <- function (copula, u, v) {
  zero <- .Machine$double.eps
  one <- 1 - .Machine$double.neg.eps
  
  u[u < zero] <- zero
  u[u > one] <- one
  v[v < zero] <- zero
  v[v > one] <- one

  rho <- copula@parameters[1]
  rho[rho == -1] <- -1 + .Machine$double.eps
  rho[rho == 1] <- 1 - .Machine$double.neg.eps
  df <- if (copula@df.fixed) copula@df else copula@parameters[2]

  r <- pt(qt(u, df+1) * sqrt(((df + qt(v, df)^2) * (1 - rho^2)) / (df+1)) + rho*qt(v, df), df)
  
  r[u <= zero | r < zero] <- zero
  r[u >= one | r > one] <- one
  
  r
}

setMethod("hinverse", "tCopula", hinverseTCopula)


hinverseClaytonCopula <- function (copula, u, v) {
  theta <- copula@parameters
  
  if (theta < .Machine$double.eps) {
    u
  } else {  
    zero <- .Machine$double.eps^0.15
    one <- 1 - .Machine$double.neg.eps^0.15
    
    u[u < zero] <- zero
    u[u > one] <- one
    v[v < zero] <- zero
    v[v > one] <- one

    r <- ((u * v^(theta+1))^(-theta/(theta+1)) + 1 - v^(-theta))^(-1/theta)
    
    r[u <= zero | r < zero] <- zero
    r[u >= one | r > one] <- one
    
    r
  }
}

setMethod("hinverse", "claytonCopula", hinverseClaytonCopula)
