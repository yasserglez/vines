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

setClass("fitVineML",
    contains = "fitVine",
    representation = representation(
        optimMethod = "character",
        optimConv = "numeric",
        startParams = "numeric",
        finalParams = "numeric",
        finalLogLik = "numeric"),
    prototype = prototype(
        method = "ml"))


logLikVine <- function (vine, data) {
  L <- function (vine, j, i, x, y) {
    # Function called by iterVine to evaluate the log-likelihood of each copula.
    copula <- vine@copulas[[j, i]]
    loglikCopula(copula@parameters, cbind(x, y), copula)
  }
  iterResult <- iterVine(vine, data, eval = L)
  sum(unlist(iterResult$evals))
}


fitVineML <- function (type, data, trees = ncol(data) - 1,
    selectCopula = function (j, i, x, y) indepCopula(),
    optimMethod = "Nelder-Mead", optimControl = list()) {
  # Compute starting values for the parameters of the copulas in the pair-copula 
  # construction following the estimation procedure described in Section 7 of 
  # Aas, K., Czado, C., Frigessi, A. and Bakken, H. Pair-copula constructions 
  # of multiple dependence. Insurance Mathematics and Economics, 2009, Vol. 44, 
  # pp. 182-198.
  selectCopulaWrapper <- function (vine, j, i, x, y) selectCopula(j, i, x, y)
  vine <- new(type, dimension = ncol(data), trees = trees,
      copulas = matrix(list(), ncol(data) - 1, ncol(data) - 1))
  vine <- iterVine(vine, data, fit = selectCopulaWrapper)$vine
  startingParams <- parameters(vine)

  if (nzchar(optimMethod) && length(startingParams) > 0) {
    # Execute the optimization method.
    lowerParams <- unlist(lapply(vine@copulas,
          function (x) if (is.null(x)) numeric(0) else x@param.lowbnd))
    upperParams <- unlist(lapply(vine@copulas, 
          function (x) if (is.null(x)) numeric(0) else x@param.upbnd))

    if (identical(optimMethod, "L-BFGS-B")) {
      lower <- lowerParams
      upper <- upperParams
    } else {
      lower <- -Inf
      upper <- Inf
    }

    L <- function (x, vine, data, lowerParams, upperParams) {
      if (all(is.finite(x) & x >= lowerParams & x <= upperParams)) {
        parameters(vine) <- x
        logLikVine(vine, data)
      } else {
        NA
      }
    }

    optimControl <- c(optimControl, fnscale = -1)
    optimResult <- optim(startingParams, L, lower = lower, upper = upper,
        method = optimMethod, control = optimControl, vine = vine, data = data, 
        lowerParams = lowerParams, upperParams = upperParams)

    parameters(vine) <- optimResult$par

    fit <- new("fitVineML", vine = vine,
        optimMethod = optimMethod,
        optimConv = optimResult$convergence,
        startParams = startingParams,
        finalParams = optimResult$par,
        finalLogLik = optimResult$value)
  } else {
    # Without parameters or optimization disabled, optimization not executed.
    fit <- new("fitVineML", vine = vine,
        optimMethod = optimMethod,
        optimConv = 0,
        startParams = startingParams,
        finalParams = startingParams,
        finalLogLik = logLikVine(vine, data))
  }

  fit
}
