# vines: R package for multivariate dependence modeling with vines
# Copyright (C) 2010-2011 Yasser González-Fernández
# Copyright (C) 2010-2011 Marta Soto
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
  eps <- .Machine$double.eps^0.5
  
  r0 <- u <= eps
  r1 <- abs(1 - u) <= eps
  s <- r0 | r1
  u <- pmax(pmin(u, 1 - eps), eps)
  v <- pmax(pmin(v, 1 - eps), eps)  
  f <- function (x, u, v) abs(h(copula, x, v) - u)
  z <- function (i) optimize(f, c(eps, 1 - eps), tol = 0.01, u = u[i], v = v[i])$minimum

  r <- sapply(seq(along = u), function(i) if (s[i]) NA else z(i))
  ifelse(r0, eps, ifelse(r1, 1 - eps, r))
}

setMethod("hinverse", "copula", hinverseCopula)


hinverseIndepCopula <- function (copula, u, v) {
  return(u)
}

setMethod("hinverse", "indepCopula", hinverseIndepCopula)


# See Aas, K., Czado, C., Frigessi, A. & Bakken, H. Pair-copula constructions 
# of multiple dependence. Insurance Mathematics and Economics, 2009, Vol. 44, 
# pp. 182-198 for the expression for the Gaussian, Student's t, Clayton and 
# Gumbel copulas.

hinverseNormalCopula <- function (copula, u, v) {
  eps <- .Machine$double.eps^0.5
  
  rho <- copula@parameters
  r0 <- u <= eps
  r1 <- abs(1 - u) <= eps
  v <- pmax(pmin(v, 1 - eps), eps)

  r <- pnorm(qnorm(u) * sqrt(1 - rho^2) + rho*qnorm(v))
  ifelse(r0, eps, ifelse(r1, 1 - eps, pmax(pmin(r, 1 - eps), eps)))
}

setMethod("hinverse", "normalCopula", hinverseNormalCopula)


hinversetCopula <- function (copula, u, v) {
  eps <- .Machine$double.eps^0.15
  
  rho <- copula@parameters
  df <- if (copula@df.fixed) copula@df else copula@parameters[2]
  r0 <- u <= eps
  r1 <- abs(1 - u) <= eps
  v <- pmax(pmin(v, 1 - eps), eps)

  r <- pt(qt(u, df+1) * sqrt(((df + qt(v, df)^2) * (1 - rho^2)) / 
              (df+1)) + rho*qt(v, df), df)
  ifelse(r0, eps, ifelse(r1, 1 - eps, pmax(pmin(r, 1 - eps), eps)))
}

setMethod("hinverse", "tCopula", hinversetCopula)


hinverseClaytonCopula <- function (copula, u, v) {
  eps <- .Machine$double.eps^0.15

  theta <- min(copula@parameters, 100)
  if (theta <= eps) return(u)
  r0 <- u <= eps
  r1 <- abs(1 - u) <= eps
  u <- pmax(pmin(u, 1 - eps), eps)
  v <- pmax(pmin(v, 1 - eps), eps)

  r <- ((u * v^(theta+1))^(-theta/(theta+1)) + 1 - v^(-theta))^(-1/theta)
  ifelse(r0, eps, ifelse(r1, 1 - eps, pmax(pmin(r, 1 - eps), eps)))
}

setMethod("hinverse", "claytonCopula", hinverseClaytonCopula)
