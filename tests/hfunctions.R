# Tests for the h-functions.

library(vines)

nvalues <- 10 # Number of values of each variable.
nparams <- 11 # Number of values of each parameter.
tolerance <- 0.01 # Tolerance checking for equality.

zeroes <- rep(.Machine$double.eps, nvalues)
ones <- rep(1 - .Machine$double.neg.eps, nvalues)
xvalues <- seq(from = zeroes[1], to = ones[1], length = nvalues)
vvalues <- seq(from = zeroes[1], to = ones[1], length = nvalues)
hargs <- merge(xvalues, vvalues)

copulas <- c(
    lapply(seq(from = -1, to = 1, length = nparams), 
        function (rho) normalCopula(rho)),
    apply(merge(seq(from = -1, to = 1, length = nparams),
            seq(from = 1, to = 30, length = nparams)), 1,
        function (p) tCopula(p[1], df = p[2], df.fixed = TRUE)),
    lapply(seq(from = .Machine$double.xmin, to = 100, length = nparams),
        function (theta) claytonCopula(theta)),
    lapply(seq(from = 1, to = 100, length = nparams),
        function (theta) gumbelCopula(theta)),    
    list(indepCopula()))

for (copula in copulas) {
  # Checking the h-function.
  
  # Check for finite return values in (0,1).
  u <- h(copula, hargs[ , 1], hargs[ , 2])
  stopifnot(all(is.finite(u)))
  stopifnot(all(u > 0 & u < 1))  

  # Check h(0, v) == 0.
  u <- h(copula, zeroes, vvalues)
  stopifnot(isTRUE(all.equal(u, zeroes, tolerance)))

  # Check h(1, v) == 1.
  u <- h(copula, ones, vvalues)
  stopifnot(isTRUE(all.equal(u, ones, tolerance)))

  # Checking the inverse of the h-functions.

  # Check for finite return values in (0,1).
  x <- hinverse(copula, u, hargs[ , 2])
  stopifnot(all(is.finite(x)))
  stopifnot(all(x > 0 & x < 1))
}
