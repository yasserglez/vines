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
  sum(unlist(iterResult$evals))
}


# Used to modify the copula before the gofCopula function is executed.
preGofCopula <- function (copula, x, y) {
  if (is(copula, "tCopula")) {
    fixedBefore <- copula@df.fixed 
    if (!fixedBefore) {
      # The gofCopula method requires the degrees-of-freedom parameter to be fixed, 
      # so it is estimated here and fixed before running the goodness-of-fit test. 
      # rho is calculated via Kendall's tau and the degrees-of-freedom is 
      # estimated by maximum likelihood with rho fixed. This follows the 
      # procedure sugested in S. Demarta and A. McNeil (2005). The t copula and 
      # related copulas. International Statistical Review 73, 111-129.
      rho <- calibKendallsTau(copula, cor(x, y, method = "kendall"))
      L <- function (df) loglikCopula(c(rho, df), cbind(x, y), copula)
      # If BFSG gives an error, try Nelder-Mead. Sometimes BFSG gives a 
      # "non-finite finite differences" error but Nelder-Mead runs OK.
      for (method in c("BFGS", "Nelder-Mead")) {
        expr <- quote(optim(copula@parameters[2], L, method = method,
                control = list(fnscale = -1)))
        optimResult <- try(suppressWarnings(eval(expr)), silent = TRUE)
        if (!inherits(optimResult, "try-error")) break
      }
      if (inherits(optimResult, "try-error")) stop(optimResult)
      copula <- tCopula(rho, df = optimResult$par, df.fixed = TRUE)
    }
    attr(copula, "df.fixed.before") <- fixedBefore
  } else if (is(copula, "claytonCopula")) {
    # Lower bound of the bivariate Clayton copula in the copula package differs
    # with Appendix B.3 of Aas, K., Czado, C., Frigessi, A. and Bakken, H. 
    # Pair-copula constructions of multiple dependence. Insurance Mathematics 
    # and Economics, 2007, Vol. 44, pp. 182-198.
    copula@param.lowbnd <- 0
  }

  copula
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
  }

  copula
}


fitVineML <- function (type, data, trees = ncol(data) - 1, copulas = list(),
    corTestMethod = "kendall", corTestSigLevel = 0.1,
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
    if (nzchar(corTestMethod)) {
      testResult <- cor.test(x, y, method = corTestMethod)
      if (testResult$p.value > corTestSigLevel) {
        selectedCopula <- indepCopula()
      }
    }

    if (is.null(selectedCopula)) {
      if (length(copulas) == 1) {
        # Only one candidate copula? Goodness-of-fit not needed.
        copula <- preGofCopula(copulas[[1]], x, y)
        # If BFSG gives an error, try Nelder-Mead. Sometimes BFSG gives a 
        # "non-finite finite differences" error but Nelder-Mead runs OK.
        for (method in c("BFGS", "Nelder-Mead")) {         
          expr <- quote(fitCopula(copula, cbind(x, y), 
                  method = gofCopulaMethod, optim.method = method,
                  estimate.variance = FALSE))
          fitResult <- try(suppressWarnings(eval(expr)), silent = TRUE)
          if (!inherits(fitResult, "try-error")) break
        }
        if (inherits(fitResult, "try-error")) stop(fitResult)
        selectedCopula <- new(class(copula), copula,
            parameters = fitResult@estimate)
        selectedCopula <- postGofCopula(selectedCopula, x, y)
      } else if (length(copulas) > 1) {
        # Goodness-of-fit tests to select the copula that better fits the data.
        pvalue <- -Inf
        for (copula in copulas) {
          copula <- preGofCopula(copula, x, y)
          # If BFSG gives an error, try Nelder-Mead. Sometimes BFSG gives a 
          # "non-finite finite differences" error but Nelder-Mead runs OK.
          for (method in c("BFGS", "Nelder-Mead")) {
            expr <- quote(gofCopula(copula, cbind(x, y), N = gofCopulaIters,
                    method = gofCopulaMethod, simulation = gofCopulaSimul,
                    optim.method = method))
            gofResult <- try(suppressWarnings(eval(expr)), silent = TRUE)
            if (!inherits(gofResult, "try-error")) break
          }
          if (inherits(gofResult, "try-error")) stop(gofResult)
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
    
    selectedCopula
  }

  # Starting values for the parameters of the copulas in the PCC following the 
  # sequential estimation procedure described in section 7 of Aas, K., 
  # Czado, C., Frigessi, A. and Bakken, H. Pair-copula constructions of multiple
  # dependence. Insurance Mathematics and Economics, 2007, Vol. 44, pp. 182-198.
  vine <- new(type, dimension = ncol(data), trees = trees,
      copulas = matrix(list(), ncol(data) - 1, ncol(data) - 1))
  iterResult <- iterVine(vine, data, fit = selectCopula)
  vine <- iterResult$vine
  startParams <- parameters(vine)
  
  if (nzchar(optimMethod) && length(startParams) > 0) {
    # Execute the optimization method.
    lowerParams <- unlist(lapply(vine@copulas,
          function (x) if (is.null(x)) numeric(0) else x@param.lowbnd))
    upperParams <- unlist(lapply(vine@copulas, 
          function (x) if (is.null(x)) numeric(0) else x@param.upbnd))

    if (optimMethod == "L-BFGS-B") {
      lower <- lowerParams
      upper <- upperParams
    } else {
      eps <- .Machine$double.eps^0.5
      lower <- -Inf + eps
      upper <- Inf - eps
    }
    
    L <- function (x, vine, data, lowerParams, upperParams) {
      if (all(x >= lowerParams) && all(x <= upperParams)) {
        parameters(vine) <- x
        logLikVine(vine, data)
      } else {
        return(NA)
      }
    }

    optimControl <- c(optimControl, fnscale = -1)
    optimResult <- optim(startParams, L, lower = lower, upper = upper,
        method = optimMethod, control = optimControl, vine = vine, data = data, 
        lowerParams = lowerParams, upperParams = upperParams)

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

  fit
}
