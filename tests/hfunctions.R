# Tests for the h functions and its inverses.

library(vines)

nValues <- 10 # Number of values of each variable.
nParams <- 11 # Number of values of each parameter.
tolerance <- 0.01 # Tolerance checking for equality.

zeroes <- rep(.Machine$double.eps, nValues)
ones <- rep(1 - .Machine$double.neg.eps, nValues)
xValues <- seq(from = zeroes[1], to = ones[1], length = nValues)
vValues <- seq(from = zeroes[1], to = ones[1], length = nValues)
hArgs <- merge(xValues, vValues)

copulas <- c(
    lapply(seq(from = -1, to = 1, length = nParams), 
        function (rho) normalCopula(rho)),
    apply(merge(seq(from = -1, to = 1, length = nParams),
            seq(from = 1, to = 500, length = nParams)), 1,
        function (p) tCopula(p[1], df = p[2], df.fixed = TRUE)),
    lapply(seq(from = 0, to = 100, length = nParams),
        function (theta) claytonCopula(theta)),
    lapply(seq(from = 1, to = 100, length = nParams),
        function (theta) gumbelCopula(theta)),    
    list(indepCopula()))

for (copula in copulas) {
  # Checking the h-function.

  # Check h(0, v) == 0.
  u <- h(copula, zeroes, vValues)
  stopifnot(isTRUE(all.equal(u, zeroes, tolerance)))

  # Check h(1, v) == 1.
  u <- h(copula, ones, vValues)
  stopifnot(isTRUE(all.equal(u, ones, tolerance)))

  # Check for finite return values in (0,1).
  u <- h(copula, hArgs[ , 1], hArgs[ , 2])
  stopifnot(all(is.finite(u)))
  stopifnot(all(u > 0 & u < 1))

  # Checking the inverse of the h-functions.

  # Check for finite return values in (0,1).
  x <- hinverse(copula, u, hArgs[ , 2])
  stopifnot(all(is.finite(x)))
  stopifnot(all(x > 0 & x < 1))

  # Cross validation.
  # stopifnot(all.equal(x, hArgs[ , 1], tolerance))
}
