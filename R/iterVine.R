# vines: R package for multivariate dependence modeling with vines
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

# Internal function that computes observations for the bivariate dependencies 
# modeled by the copulas in the vine (using the observations of the original 
# variables in the data matrix) and executes the f function with the vine, the 
# indexes of the copula and the observations of each variable as arguments.

setGeneric("iterVine", 
    function (vine, data, fit = NULL, eval = NULL) {
      if (identical(vine@trees, 0)) {
        # Vine without trees, nothing to iterate for.
        list(vine = vine, evals = list())
      } else {
        standardGeneric("iterVine")
      }
    },
    signature = "vine")


iterCVine <- function (vine, data, fit = NULL, eval = NULL) {
  # The implementation of this function is based on the Algorithm 3 described in 
  # Aas, K., Czado, C., Frigessi, A. & Bakken, H. Pair-copula constructions of 
  # multiple dependence. Insurance Mathematics and Economics, 2009, Vol. 44, 
  # pp. 182-198.
  
  # The indexes of the second dimention of the v array differs with the 
  # indexes of the first dimention of the v array in Algorithm 3 because of 
  # GNU R 1-based indexing.
  
  # This implementation avoids evaluating the h-functions beyond the last tree 
  # of the vine that represents dependence (given by the trees slot of the vine) 
  # because the h-functions of the Independence copula always return the value 
  # of its first argument.

  evals <- list()
  d <- vine@dimension
  v <- array(NA, c(nrow(data), d - 1, d))

  for (i in seq(length = d)) {
    v[ , 1, i] <- data[ , i]
  }
  for (j in seq(length = vine@trees)) {
    for (i in seq(length = d - j)) {
      if (!is.null(fit)) {
        vine@copulas[[j, i]] <- fit(vine, j, i, v[ , j, 1], v[ , j, i+1])
      }
      if (!is.null(eval)) {
        evals <- c(evals, list(eval(vine, j, i, v[ , j, 1], v[ , j, i+1])))
      }
    }

    if (identical(j, vine@trees)) break

    # Compute observations for the next tree.
    for (i in seq(length = d - j)) {
      v[ , j+1, i] <- h(vine@copulas[[j, i]], v[ , j, i+1], v[ , j, 1])
    }
  }

  list(vine = vine, evals = evals)
}

setMethod("iterVine", "CVine", iterCVine)


iterDVine <- function (vine, data, fit = NULL, eval = NULL) {
  # The implementation of this function is based on the Algorithm 4 described in 
  # Aas, K., Czado, C., Frigessi, A. & Bakken, H. Pair-copula constructions of 
  # multiple dependence. Insurance Mathematics and Economics, 2009, Vol. 44, 
  # pp. 182-198.
  
  # The indexes of the second dimention of the v array differs with the 
  # indexes of the first dimention of the v array in Algorithm 4 because of 
  # GNU R 1-based indexing.
  
  # This implementation avoids evaluating the h-functions beyond the last tree 
  # of the vine that represents dependence (given by the trees slot of the vine) 
  # because the h-functions of the Independence copula always return the value 
  # of its first argument.

  evals <- list()
  d <- vine@dimension
  v <- array(NA, c(nrow(data), d, max(2 * d - 4, d)))

  for (i in seq(length = d)) {
    v[ , 1, i] <- data[ , i]
  }
  for (i in seq(length = d - 1)) {
    if (!is.null(fit)) {
      vine@copulas[[1, i]] <- fit(vine, 1, i, v[ , 1, i], v[ , 1, i+1])
    }
    if (!is.null(eval)) {
      evals <- c(evals, list(eval(vine, 1, i, v[ , 1, i], v[ , 1, i+1])))    
    }
  }
  v[ , 2, 1] <- h(vine@copulas[[1, 1]], v[ , 1, 1], v[ , 1, 2])
  for (k in seq(length = max(d - 3, 0))) {
    v[ , 2, 2*k] <- h(vine@copulas[[1, k+1]], v[ , 1, k+2], v[ , 1, k+1])
    v[ , 2, 2*k+1] <- h(vine@copulas[[1, k+1]], v[ , 1, k+1], v[ , 1, k+2])
  }
  v[ , 2, 2*d-4] <- h(vine@copulas[[1, d-1]], v[ , 1, d], v[ , 1, d-1])
  for (j in seq(from = 2, length = vine@trees - 1)) {
    for (i in seq(length = d - j)) {
      if (!is.null(fit)) {
        vine@copulas[[j, i]] <- fit(vine, j, i, v[ , j, 2*i-1], v[ , j, 2*i])
      }
      if (!is.null(eval)) {
        evals <- c(evals, list(eval(vine, j, i, v[ , j, 2*i-1], v[ , j, 2*i])))
      }
    }

    if (identical(j, vine@trees)) break

    # Compute observations for the next tree.
    v[ , j+1, 1] <- h(vine@copulas[[j, 1]], v[ , j, 1], v[ , j, 2])
    if (d > 4) {
      for (i in seq(length = d - j - 2)) {
        v[ , j+1, 2*i] <- h(vine@copulas[[j, i+1]], v[ , j, 2*i+2], v[ , j, 2*i+1])
        v[ , j+1, 2*i+1] <- h(vine@copulas[[j, i+1]], v[ , j, 2*i+1], v[ , j, 2*i+2])
      }
    }
    v[ , j+1, 2*d-2*j-2] <- h(vine@copulas[[j, d-j]], v[ , j, 2*d-2*j], v[ , j, 2*d-2*j-1])
  }

  list(vine = vine, evals = evals)
}

setMethod("iterVine", "DVine", iterDVine)
