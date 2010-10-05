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

orderingVineGreedy <- function (type, data, according = "kendall") {
  if (according %in% c("pearson", "kendall", "spearman")) {
    values <- abs(cor(data, method = according))
    diag(values) <- -Inf
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
  
  if (type == "DVine") {
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
  } else if (type == "CVine") {
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


orderingVine <- function (type, data, method = "greedy", ...) {
  if (type %in% c("CVine", "DVine") && method == "greedy") {
    orderingVineGreedy(type, data, ...)
  } else {
    stop(paste("invalid", sQuote(method), "ordering method for", type))
  }
}

orderingCVine <- function (data, method = "greedy", ...) {
  orderingVine("CVine", data, method, ...)
}

orderingDVine <- function (data, method = "greedy", ...) {
  orderingVine("DVine", data, method, ...)
}
