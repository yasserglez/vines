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

parameters <- function (vine) {
    f <- function (x) if (is.null(x)) numeric(0) else x@parameters
    unlist(lapply(vine@copulas, f))
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
    vine
}
