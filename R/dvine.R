# vines: Multivariate Dependence Modeling with Vines
# Copyright (C) 2010-2013 Yasser Gonzalez-Fernandez <ygonzalezfernandez@gmail.com>
# Copyright (C) 2010-2013 Marta Soto <mrosa@icimaf.cu>
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

setGeneric("dvine",
        function (vine, u) standardGeneric("dvine"),
        signature = "vine")


dCVineDVine <- function (vine, u) {
    if (is.vector(u)) u <- matrix(u, nrow = 1)

    if (vine@trees == 0) {
        rep(1, nrow(u)) # The product of the uniform marginal densities.
    } else {
        evalCopula <- function (vine, j, i, x, y) {
            dCopula(cbind(x, y), vine@copulas[[j, i]], log = TRUE)
        }
        iterResult <- vineIter(vine, u, evalCopula = evalCopula)
        exp(apply(matrix(unlist(iterResult$evals), nrow(u)), 1, sum))
    }
}

setMethod("dvine", "CVine", dCVineDVine)
setMethod("dvine", "DVine", dCVineDVine)
