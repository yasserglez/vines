# vines: R package for multivariate dependence modeling with vines
# Copyright (C) 2010 Yasser González-Fernández <ygonzalezfernandez@gmail.com>
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

setClass("fitVine",
    representation = representation(
        vine = "Vine",
        observations = "numeric",
        method = "character"))


showFitVine <- function (object) {
  cat("Vine Inference\n\n")
  cat("Method:", object@method, "\n")
  cat("Vine type:", object@vine@type, "\n")
  cat("Dimension:", object@vine@dimension, "\n")
  cat("Observations:", object@observations, "\n")
}

setMethod("show", "fitVine", showFitVine)


fitVine <- function (type, data, method = "ml", ...) {
  if (type %in% c("CVine", "DVine") && identical(method, "ml")) {
    fitVineML(type, data, ...)
  } else {
    stop("invalid fit method ", dQuote(method), " for ", dQuote(type))
  }
}

fitCVine <- function (data, method = "ml", ...) {
  fitVine("CVine", data, method, ...)
}

fitDVine <- function (data, method = "ml", ...) {
  fitVine("DVine", data, method, ...)
}
