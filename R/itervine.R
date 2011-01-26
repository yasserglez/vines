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

setGeneric("iterVine",
        function (vine, data, fit = NULL, eval = NULL)
            standardGeneric("iterVine"),
        signature = "vine")


iterCVine <- function (vine, data, fit = NULL, eval = NULL) {
    .Call(C_iterCVine, vine, data, fit, eval)
}

setMethod("iterVine", "CVine", iterCVine)


iterDVine <- function (vine, data, fit = NULL, eval = NULL) {
    .Call(C_iterDVine, vine, data, fit, eval)
}

setMethod("iterVine", "DVine", iterDVine)
