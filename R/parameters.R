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

# Functions to get and set the parameters of all the copulas in a vine. 
# The parameters function returns a numeric vector with the parameters of all 
# the copulas concatenated. The replacement function parameters<- receives a 
# single vector with the parameters of all the copulas concatenated and updates 
# the copulas matrix of the vine. The parameters of the copulas are 
# concatenated iterating the copulas matrix by columns.

parameters <- function (vine) {
  f <- function (x) if (is.null(x)) numeric(0) else x@parameters
  return(unlist(lapply(vine@copulas, f)))
}


`parameters<-` <- function (vine, value) {
  k <- 1
  for (i in seq(along = vine@copulas)) {
    if (!is.null(vine@copulas[[i]])) {
      n <- length(vine@copulas[[i]]@parameters)
      if (n > 0) {
        params <- value[seq(from = k, to = k + n - 1)]
        vine@copulas[[i]]@parameters <- params
        k <- k + n
      }
    }
  }
  return(vine)
}
