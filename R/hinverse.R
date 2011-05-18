# vines: Multivariate Dependence Modeling with Vines
# Copyright (C) 2010, 2011 Yasser González-Fernández <ygf@icmf.inf.cu>
# Copyright (C) 2010, 2011 Marta Soto <mrosa@icmf.inf.cu>
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
    skip <- r0 | r1
    u <- pmax(pmin(u, 1-eps), eps)
    v <- pmax(pmin(v, 1-eps), eps)
    f <- function (x, u, v, copula) h(copula, x, v) - u
    r <- sapply(seq(along = u),
            function (i) {
                if (skip[i]) {
                    NA
                } else {
                    uniroot(f, lower = eps, upper = 1-eps,
                            f.lower = -u[i], f.upper = 1-u[i], tol = 0.01,
                            copula = copula, u = u[i], v = v[i])$root
                }
            })
    ifelse(r0, eps, ifelse(r1, 1 - eps, r))
}

setMethod("hinverse", "copula", hinverseCopula)


hinverseIndepCopula <- function (copula, u, v) {
    .Call(C_hinverseIndepCopula, u, v)
}

setMethod("hinverse", "indepCopula", hinverseIndepCopula)


hinverseNormalCopula <- function (copula, u, v) {
    rho <- copula@parameters
    .Call(C_hinverseNormalCopula, rho, u, v)
}

setMethod("hinverse", "normalCopula", hinverseNormalCopula)


hinverseTCopula <- function (copula, u, v) {
    rho <- copula@parameters
    df <- if (copula@df.fixed) copula@df else copula@parameters[2]
    .Call(C_hinverseTCopula, rho, df, u, v)
}

setMethod("hinverse", "tCopula", hinverseTCopula)


hinverseClaytonCopula <- function (copula, u, v) {
    theta <- min(copula@parameters, 100)
    .Call(C_hinverseClaytonCopula, theta, u, v)
}

setMethod("hinverse", "claytonCopula", hinverseClaytonCopula)


hinverseFrankCopula <- function (copula, u, v) {
    theta <- max(min(copula@parameters, 100), -100)
    .Call(C_hinverseFrankCopula, theta, u, v)
}

setMethod("hinverse", "frankCopula", hinverseFrankCopula)
