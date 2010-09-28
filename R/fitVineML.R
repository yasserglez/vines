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
  # Function called by iterVine to evaluate the log-likelihood of each copula.
  L <- function (vine, i, j, x, y) {
    copula <- vine@copulas[[i, j]]
    loglikCopula(copula@parameters, cbind(x, y), copula)
  }
  iterResult <- iterVine(vine, data, eval = L)
  return(sum(unlist(iterResult$evals)))
}


# Used to modify the copula before the gofCopula function is executed.
preGofCopula <- function (copula, x, y) {
  if (is(copula, "tCopula")) {
    fixedBefore <- copula@df.fixed 
    if (!fixedBefore) {
      # The gofCopula method requires the dof parameter to be fixed. We 
      # estimate and fix this parameters before running the GoF test.
      fitResult <- fitCopula(tCopula(0), cbind(x, y), method = "ml", 
          start = c(0, 4), estimate.variance = FALSE)
      rho <- fitResult@copula@parameters[1]
      df <- fitResult@copula@parameters[2]
      copula <- tCopula(rho, df = df, df.fixed = TRUE)
    }
    attr(copula, "df.fixed.before") <- fixedBefore
    copula
  } else {
    # Return the copula without modifications.
    copula
  }
}


# Used to modify the copula after the gofCopula function is executed.
postGofCopula <- function (copula, x, y) {
  if (is(copula, "tCopula")) {
    fixedBefore <- attr(copula, "df.fixed.before")
    if (!fixedBefore) {
      # Unfix the degrees-of-freedom so it is considered as a variable
      # during the optimization of the log-likelihood of the vine.
      rho <- copula@parameters[1]
      df <- copula@df
      copula <- tCopula(rho, df = df, df.fixed = FALSE)
    }
    attr(copula, "df.fixed.before") <- NULL
    copula
  } else {
    # Return the copula without modifications.
    copula
  }
}


fitVineML <- function (type, data, trees = ncol(data) - 1, copulas = list(), 
    corTestMethod = "kendall", corTestSigLevel = 0.05, 
    gofCopulaIters = 1000, gofCopulaMethod = "mpl", gofCopulaSimul = "mult", 
    optimMethod = "Nelder-Mead", optimControl = list()) {
  
  if (is.null(corTestMethod)) corTestMethod <- ""
  if (is.null(optimMethod)) optimMethod <- ""
  
  # Function called by iterVine to select the bivariate copula that better 
  # fits the bivariate data from the given list of copulas.
  selectCopula <- function (vine, i, j, x, y) {
    selectedCopula <- NULL
    # Test for association. The Independence copula will be used if there
    # is not enough evidence of association between the variables.
    if (corTestMethod != "") {
      testResult <- cor.test(x, y, method = corTestMethod)
      if (testResult$p.value > corTestSigLevel) {
        selectedCopula <- indepCopula()
      }
    }

    if (is.null(selectedCopula)) {
      if (length(copulas) == 1) {
        # Only one candidate copula? Goodness-of-fit not needed.
        copula <- preGofCopula(copulas[[1]], x, y)
        fitResult <- fitCopula(copula, cbind(x, y), method = gofCopulaMethod,
            optim.method = "Nelder-Mead", estimate.variance = FALSE)
        selectedCopula <- new(class(copula), copula,
            parameters = fitResult@estimate)
        selectedCopula <- postGofCopula(selectedCopula, x, y)
      } else if (length(copulas) > 1) {
        # Goodness-of-fit tests to select the copula that better fits the data.
        pvalue <- -Inf
        for (copula in copulas) {
          copula <- preGofCopula(copula, x, y)
          gofResult <- gofCopula(copula, cbind(x, y), N = gofCopulaIters,
              optim.method = "Nelder-Mead", method = gofCopulaMethod, 
              simulation = gofCopulaSimul)
          if (gofResult$pvalue > pvalue) {
            pvalue <- gofResult$pvalue
            selectedCopula <- new(class(copula), copula, 
                parameters = gofResult$parameters)
            selectedCopula <- postGofCopula(selectedCopula, x, y)
          }
        }
      }
    }

    # Without candidate copulas? Independence vine. 
    if (is.null(selectedCopula)) {
      selectedCopula <- indepCopula()
    }
    
    return(selectedCopula)
  }

  # Starting values for the parameters of the copulas in the PCC following the 
  # sequential estimation procedure described in section 7 of Aas, K., 
  # Czado, C., Frigessi, A. & Bakken, H. Pair-copula constructions of multiple
  # dependence. Insurance Mathematics and Economics, 2007, Vol. 44, pp. 182-198.
  vine <- new(type, dimension = ncol(data), trees = trees,
      copulas = matrix(list(), ncol(data), ncol(data)))
  iterResult <- iterVine(vine, data, fit = selectCopula)
  vine <- iterResult$vine
  startParams <- parameters(vine)
  
  if (optimMethod != "" && length(startParams) > 0) {
    # Execute the optimization method.
    if (optimMethod == "L-BFGS-B") {
      lower <- unlist(lapply(vine@copulas, 
            function (x) if (is.null(x)) numeric(0) else 
                x@param.lowbnd + (.Machine$double.eps ^ 0.15)))
      upper <- unlist(lapply(vine@copulas, 
            function (x) if (is.null(x)) numeric(0) else 
                x@param.upbnd - (.Machine$double.eps ^ 0.15)))
    } else {
      lower <- -Inf
      upper <- Inf
    }

    L <- function (x, vine, data) {
      parameters(vine) <- x
      logLikVine(vine, data) 
    }
    optimControl <- c(optimControl, fnscale = -1)
    optimResult <- optim(startParams, L, lower = lower, upper = upper,
        method = optimMethod, control = optimControl,
        vine = vine, data = data)
    
    parameters(vine) <- optimResult$par

    fit <- new("fitVineML", vine = vine,
        optimMethod = optimMethod,
        optimConv = optimResult$convergence,
        startParams = startParams,
        finalParams = optimResult$par,
        finalLogLik = optimResult$value)
  } else {
    # Without parameters or optimization disabled, optimization not executed.
    fit <- new("fitVineML", vine = vine,
        optimMethod = optimMethod,
        optimConv = 0,
        startParams = startParams,
        finalParams = startParams,
        finalLogLik = logLikVine(vine, data))
  }

  return(fit)
}
