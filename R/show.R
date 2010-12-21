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

setMethod("toString", "copula",
    function (x, ...) {
      if (length(x@parameters) > 0) {
        parts <- character(0)
        for (k in seq(along = x@parameters)) {
          parts <- c(parts, paste(x@param.names[k], "=", x@parameters[k]))
        }
        parts <- paste(parts, collapse = ", ")
        parameters <- paste(" (", parts, ")", sep = "")
      } else {
        parameters <- ""
      }
      paste(sub('[[:space:]]+$', '', x@message), parameters, sep = "")
    })


showVine <- function (object) {
  cat("Vine\n\n")
  cat("Type:", object@type, "\n")
  cat("Dimension:", object@dimension, "\n")
  cat("Dependency trees:", object@trees, "\n")
}

setMethod("show", "Vine", showVine)


showCVine <- function (object) {
  showVine(object)
  if (length(object@dimensionNames) > 0) {
    dimNames <- object@dimensionNames
  } else {
    dimNames <- as.character(seq(length = object@dimension))
  }
  for (j in seq(length = object@trees)) {
    copulaLabels <- character(0)
    copulaDescs <- character(0)

    for (i in seq(length = object@dimension - j)) {
      conditioned <- paste(dimNames[j], dimNames[j + i], sep = ",")
      conditioning <- paste(dimNames[seq(length = j - 1)], collapse = ",")
      copulaLabel <- paste(conditioned,
          if (j > 1) paste("|", conditioning, sep = "")
              else character(0), 
          sep = "")

      copulaLabels <- c(copulaLabels, copulaLabel)
      copulaDescs <- c(copulaDescs, toString(object@copulas[[j, i]]))
    }

    cat("\n")
    copulaLabels <- format(copulaLabels)
    for (k in seq(along = copulaLabels)) {
      cat(copulaLabels[k], ": ", copulaDescs[k], "\n", sep = "")
    }

    if (j == object@trees) break
  }
}

setMethod("show", "CVine", showCVine)


showDVine <- function (object) {
  showVine(object)
  if (length(object@dimensionNames) > 0) {
    dimNames <- object@dimensionNames
  } else {
    dimNames <- as.character(seq(length = object@dimension))
  }
  for (j in seq(length = object@trees)) {
    copulaLabels <- character(0)
    copulaDescs <- character(0)

    for (i in seq(length = object@dimension - j)) {
      conditioned <- paste(dimNames[i], dimNames[i + j], sep = ",")
      conditioning <- paste(dimNames[seq(from = i + 1, to = i + j - 1)], collapse = ",")
      copulaLabel <- paste(conditioned,
              if (j > 1) paste("|", conditioning, sep = "")
              else character(0), 
          sep = "")

      copulaLabels <- c(copulaLabels, copulaLabel)
      copulaDescs <- c(copulaDescs, toString(object@copulas[[j, i]]))
    }

    cat("\n")
    copulaLabels <- format(copulaLabels)
    for (k in seq(along = copulaLabels)) {
      cat(copulaLabels[k], ": ", copulaDescs[k], "\n", sep = "")
    }

    if (j == object@trees) break
  }
}

setMethod("show", "DVine", showDVine)
