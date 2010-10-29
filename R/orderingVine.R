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

orderingVineCanonical <- function (type, data) {
  seq(from = 1, to = ncol(data))
}


orderingVineGreedy <- function (type, data, according = "kendall") {
  if (according %in% c("pearson", "kendall", "spearman")) {
    # Calculate the value of the given measure of dependence between each
    # pair of variables and couple the variables with the largest modular values.
    values <- abs(cor(data, method = according))
    diag(values) <- -Inf
  } else if (according %in% c("df")) {
    # Fit bivariate t copulas to each pair of variables and copule the variables
    # with the smaller degrees of freedom.
    values <- matrix(-Inf, ncol(data), ncol(data))
    for (i in seq(length = ncol(data))) {
      for (j in seq(length = max(i - 1, 0))) {
        x <- data[ , i]
        y <- data[ , j]
        copula <- tCopula(0)
        rho <- calibKendallsTau(copula, cor(x, y, method = "kendall"))
        L <- function (df) loglikCopula(c(rho, df), cbind(x, y), copula)
        df <- optim(copula@parameters[2], L, method = "BFGS",
            control = list(fnscale = -1))$par
        values[i, j] <- -df
        values[j, i] <- -df
      }
    }
  } else {
    stop('invalid argument of the according argument')
  }

  # Try to couple the pairs with the maximum value in the values matrix, 

  n <- ncol(data)
  k <- which.max(values)
  i <- row(values)[k]
  j <- col(values)[k]
  values[i, j] <- -Inf
  values[j, i] <- -Inf
  
  if (identical(type, "DVine")) {
    ordering <- c(i, j)
    while (length(ordering) < n) {
      i <- ordering[1]
      j <- ordering[length(ordering)]
      if (max(values[i, ]) >= max(values[j, ])) {
        ii <- which.max(values[i, ])
        for (k in ordering) {
          values[k, ii] <- -Inf
          values[ii, k] <- -Inf
        }
        ordering <- c(ii, ordering)
      } else {
        jj <- which.max(values[j, ])
        for (k in ordering) {
          values[k, jj] <- -Inf
          values[jj, k] <- -Inf
        }
        ordering <- c(ordering, jj)
      }
    }
  } else if (identical(type, "CVine")) {
    if (max(values[i, ]) >= max(values[j, ])) {
      ii <- which.max(values[i, ])
      ordering <- c(i, ii, j, seq(to = n)[c(-i, -ii, -j)])
    } else {
      jj <- which.max(values[j, ])
      ordering <- c(j, jj, i, seq(to = n)[c(-j, -jj, -i)])
    }
  }
  
  ordering
}


orderingVine <- function (type, data, method = "canonical", ...) {
  if (type %in% c("CVine", "DVine") && identical(method, "greedy")) {
    orderingVineGreedy(type, data, ...)
  } else if (type %in% c("CVine", "DVine") && identical(method, "canonical")) {
    orderingVineCanonical(type, data)
  } else {
    stop(paste("invalid", sQuote(method), "ordering method for", type))
  }
}

orderingCVine <- function (data, method = "canonical", ...) {
  orderingVine("CVine", data, method, ...)
}

orderingDVine <- function (data, method = "canonical", ...) {
  orderingVine("DVine", data, method, ...)
}
