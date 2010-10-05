# edas: GNU R package for Estimation of Distribution Algorithms
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

# The lower bound of the parameter of the bivariate Clayton copula in the copula
# package differs with Appendix B.3 of Aas, K., Czado, C., Frigessi, A. and
# Bakken, H. Pair-copula constructions of multiple dependence. Insurance 
# Mathematics and Economics, 2007, Vol. 44, pp. 182-198.

calibKendallsTauClaytonCopula <- function (copula, tau) {
  max(copula:::calibKendallsTauClaytonCopula(copula, tau), 
      copula@param.lowbnd + .Machine$double.eps^0.5)
}


calibSpearmansRhoClaytonCopula <- function (copula, rho) {
  max(copula:::calibSpearmansRhoClaytonCopula(copula, rho), 
      copula@param.lowbnd + .Machine$double.eps^0.5)
}


.onLoad <- function (libname, pkgname) {
  setMethod("calibKendallsTau", signature("claytonCopula"), 
      calibKendallsTauClaytonCopula)
  
  setMethod("calibSpearmansRho", signature("claytonCopula"),
      calibSpearmansRhoClaytonCopula)  
}
