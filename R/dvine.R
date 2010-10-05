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

setGeneric("dvine",
    function (vine, u) {
      if (is.vector(u)) u <- matrix(u, nrow = 1)
      if (vine@trees == 0) {
        # Vine without trees, calculate the product of the margins.
        apply(u, 1, prod)
      } else {
        standardGeneric("dvine")
      }
    },
    signature = "vine")


dCVineDVine <- function (vine, u) {
  # Function called by iterVine to evaluate the density of each copula.
  copulaDensity <- function (vine, i, j, x, y) {
    dcopula(vine@copulas[[i, j]], cbind(x, y))
  }
  iterResult <- iterVine(vine, u, eval = copulaDensity)
  copulasProd <- apply(matrix(unlist(iterResult$evals), nrow(u)), 1, prod)
  marginsProd <- apply(u, 1, prod)
  copulasProd * marginsProd
}

setMethod("dvine", "CVine", dCVineDVine)

setMethod("dvine", "DVine", dCVineDVine)
