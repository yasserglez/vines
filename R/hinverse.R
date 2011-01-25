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
    .Call(C_hinverseCopula, copula, u, v)
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
    eps <- .Machine$double.eps^0.15
    theta <- min(copula@parameters, 100)
    .Call(C_hinverseClaytonCopula, theta, u, v)
}

setMethod("hinverse", "claytonCopula", hinverseClaytonCopula)
