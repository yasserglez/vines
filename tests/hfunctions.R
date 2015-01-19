library("vines")

n <- 7 # Number of values of each variable.
np <- 5 # Number of values of each parameter.
tol <- 0.01 # Ignore differences smaller than tol.

X <- seq(from = 0.25, to = 0.75, length = n)
V <- seq(from = 0.25, to = 0.75, length = n)
XV <- merge(X, V)

copulas <- c(
    lapply(seq(from = -1, to = 1, length = np),
        function (rho) normalCopula(rho)),
    apply(merge(seq(from = -1, to = 1, length = np),
            seq(from = 1, to = 30, length = np)), 1,
        function (p) tCopula(p[1], df = as.integer(p[2]), df.fixed = TRUE)),
    lapply(seq(from = .Machine$double.eps^0.5, to = 3, length = np),
        function (theta) claytonCopula(theta, use.indepC = "FALSE")),
    lapply(seq(from = 1, to = 5, length = np),
        function (theta) gumbelCopula(theta, use.indepC = "FALSE")),
    lapply(seq(from = -1, to = 1, length = np),
        function (theta) fgmCopula(theta)),
    lapply(seq(from = 1, to = 5, length = np),
        function (theta) galambosCopula(theta)),
    lapply(seq(from = -10, to = 10, length = np),
        function (theta) frankCopula(theta, use.indepC = "FALSE")),
    list(indepCopula()))

for (copula in copulas) {
    # Validate the h-function.
    u <- h(copula, XV[ , 1], XV[ , 2])
    uu <- vines:::hCopula(copula, XV[ , 1], XV[ , 2])
    stopifnot(isTRUE(all.equal(u, uu, tolerance = tol)))

    # Validate the inverse of the h-function.
    x <- hinverse(copula, u, XV[ , 2])
    xx <- vines:::hinverseCopula(copula, u, XV[ , 2])
    stopifnot(isTRUE(all.equal(x, xx, tolerance = tol)))
}


X <- seq(from = 0, to = 1, length = n)
V <- seq(from = 0, to = 1, length = n)
XV <- merge(X, V)

copulas <- c(
    lapply(seq(from = -1, to = 1, length = np),
        function (rho) normalCopula(rho)),
    apply(merge(seq(from = -1, to = 1, length = np),
            seq(from = 1, to = 30, length = np)), 1,
        function (p) tCopula(p[1], df = as.integer(p[2]), df.fixed = TRUE)),
    lapply(seq(from = .Machine$double.eps^0.5, to = 100, length = np),
        function (theta) claytonCopula(theta, use.indepC = "FALSE")),
    lapply(seq(from = 1, to = 100, length = np),
        function (theta) gumbelCopula(theta, use.indepC = "FALSE")),
    lapply(seq(from = -1, to = 1, length = np),
        function (theta) fgmCopula(theta)),
    lapply(seq(from = 1, to = 25, length = np),
        function (theta) galambosCopula(theta)),
    lapply(seq(from = -100, to = 100, length = np),
        function (theta) frankCopula(theta, use.indepC = "FALSE")),
    list(indepCopula()))

for (copula in copulas) {
    # Check h(0, v) == 0.
    u <- h(copula, rep(0, n), V)
    stopifnot(isTRUE(all.equal(u, rep(0, n), tolerance = tol)))

    # Check h(1, v) == 1.
    u <- h(copula, rep(1, n), V)
    stopifnot(isTRUE(all.equal(u, rep(1, n), tolerance = tol)))

    # Check for finite return values in [0,1].
    u <- h(copula, XV[ , 1], XV[ , 2])
    stopifnot(all(is.finite(u)))
    stopifnot(all(u >= 0 & u <= 1))

    # Check for finite return values in [0,1].
    x <- hinverse(copula, u, XV[ , 2])
    stopifnot(all(is.finite(x)))
    stopifnot(all(x >= 0 & x <= 1))

    # Check hinverse(0, v) == 0.
    x <- hinverse(copula, rep(0, n), V)
    stopifnot(isTRUE(all.equal(x, rep(0, n), tolerance = tol)))

    # Check hinverse(1, v) == 1.
    x <- hinverse(copula, rep(1, n), V)
    stopifnot(isTRUE(all.equal(x, rep(1, n), tolerance = tol)))
}
