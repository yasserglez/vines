# vines: R package for multivariate dependence modeling with vines
# Copyright (C) 2010, 2011 Yasser González-Fernández
# Copyright (C) 2010, 2011 Marta Soto
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
  eps <- .Machine$double.eps^0.5
  
  envir <- new.env()
  assign("copula", copula, envir)
  assign("x", pmax(pmin(x, 1 - eps), eps), envir)
  assign("v", pmax(pmin(v, 1 - eps), eps), envir)
  d <- numericDeriv(quote(pcopula(copula, cbind(x, v))), "v", envir)
  r <- diag(attr(d, "gradient"))
  pmax(pmin(r, 1 - eps), eps)
}

setMethod("h", "copula", hCopula)


hIndepCopula <- function (copula, x, v) {
  .Call(C_hIndepCopula, x, v)
}

setMethod("h", "indepCopula", hIndepCopula)


# See Aas, K., Czado, C., Frigessi, A. & Bakken, H. Pair-copula constructions 
# of multiple dependence. Insurance Mathematics and Economics, 2009, Vol. 44, 
# pp. 182-198 for the expression for the Gaussian, Student's t, Clayton and 
# Gumbel copulas.

hNormalCopula <- function (copula, x, v) {
  rho <- copula@parameters
  .Call(C_hNormalCopula, rho, x, v)
}

setMethod("h", "normalCopula", hNormalCopula)


hTCopula <- function (copula, x, v) {
  rho <- copula@parameters
  df <- if (copula@df.fixed) copula@df else copula@parameters[2]
  .Call(C_hTCopula, rho, df, x, v)
}

setMethod("h", "tCopula", hTCopula)


hClaytonCopula <- function (copula, x, v) {
  theta <- min(copula@parameters, 100)
  .Call(C_hClaytonCopula, theta, x, v)
}

setMethod("h", "claytonCopula", hClaytonCopula)


hGumbelCopula <- function (copula, x, v) {
  theta <- min(copula@parameters, 100)
  .Call(C_hGumbelCopula, theta, x, v)
}

setMethod("h", "gumbelCopula", hGumbelCopula)
