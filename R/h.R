# vines: R package for multivariate dependence modeling with vines
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
    function (copula, x, v) standardGeneric("h"),
    signature = "copula")


hCopula <- function (copula, x, v) {
  # Last resort is to evaluate the h-function numerically.
  zero <- .Machine$double.eps
  one <- 1 - .Machine$double.neg.eps
  
  x[x < zero] <- zero
  x[x > one] <- one
  v[v < zero] <- zero
  v[v > one] <- one    

  e <- new.env()
  assign("x", x, envir = e)
  assign("v", v, envir = e)
  d <- numericDeriv(quote(pcopula(copula, cbind(x, v))), c("v"), e)
  r <- diag(attr(d, "gradient"))
  
  r[x <= zero | r < zero] <- zero
  r[x >= one | r > one] <- one

  r
}

setMethod("h", "copula", hCopula)


hIndepCopula <- function (copula, x, v) {
  x
}

setMethod("h", "indepCopula", hIndepCopula)


# See Aas, K., Czado, C., Frigessi, A. & Bakken, H. Pair-copula constructions 
# of multiple dependence. Insurance Mathematics and Economics, 2009, Vol. 44, 
# pp. 182-198 for the expression for the Gaussian, Student's t, Clayton and 
# Gumbel copulas.

hNormalCopula <- function (copula, x, v) {
  zero <- .Machine$double.eps
  one <- 1 - .Machine$double.neg.eps

  x[x < zero] <- zero
  x[x > one] <- one
  v[v < zero] <- zero
  v[v > one] <- one  

  rho <- copula@parameters
  rho[rho == -1] <- -1 + .Machine$double.eps
  rho[rho == 1] <- 1 - .Machine$double.neg.eps

  r <- pnorm((qnorm(x) - rho*qnorm(v)) / sqrt(1 - rho^2))

  r[x <= zero | r < zero] <- zero
  r[x >= one | r > one] <- one

  r  
}

setMethod("h", "normalCopula", hNormalCopula)


hTCopula <- function (copula, x, v) {
  zero <- .Machine$double.eps
  one <- 1 - .Machine$double.neg.eps

  x[x < zero] <- zero
  x[x > one] <- one
  v[v < zero] <- zero
  v[v > one] <- one    

  rho <- copula@parameters[1]
  rho[rho == -1] <- -1 + .Machine$double.eps
  rho[rho == 1] <- 1 - .Machine$double.neg.eps  
  df <- if (copula@df.fixed) copula@df else copula@parameters[2]

  r <- pt((qt(x, df) - rho*qt(v, df)) / sqrt(((df + qt(v, df)^2) * (1 - rho^2)) / (df+1)), df+1)

  r[x <= zero | r < zero] <- zero
  r[x >= one | r > one] <- one
  
  r
}

setMethod("h", "tCopula", hTCopula)


hClaytonCopula <- function (copula, x, v) {
  theta <- copula@parameters

  if (theta < .Machine$double.eps) {
    x
  } else {
    zero <- .Machine$double.eps^0.15
    one <- 1 - .Machine$double.neg.eps^0.15

    x[x < zero] <- zero
    x[x > one] <- one
    v[v < zero] <- zero
    v[v > one] <- one  
  
    r <- v^(-theta-1) * (x^(-theta) + v^(-theta) - 1)^(-1-1/theta)
  
    r[x <= zero | r < zero] <- zero
    r[x >= one | r > one] <- one
  
    r
  }
}

setMethod("h", "claytonCopula", hClaytonCopula)


hGumbelCopula <- function (copula, x, v) {
  zero <- .Machine$double.eps
  one <- 1 - .Machine$double.neg.eps
  
  x[x < zero] <- zero
  x[x > one] <- one
  v[v < zero] <- zero
  v[v > one] <- one  
  
  theta <- copula@parameters

  r <- pcopula(copula, cbind(x, v)) * 1/v * (-log(v))^(theta-1) *
      ((-log(x))^theta + (-log(v))^theta)^(1/theta-1)

  r[x <= zero | r < zero] <- zero
  r[x >= one | r > one] <- one

  r
}

setMethod("h", "gumbelCopula", hGumbelCopula)
