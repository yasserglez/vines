library("vines")

N <- 7 # Number of values of each variable.
P <- 5 # Number of values of each parameter.
T <- 0.01 # Tolerance checking for equality.

X <- seq(from = 0.25, to = 0.75, length = N)
V <- seq(from = 0.25, to = 0.75, length = N)
XV <- merge(X, V)

copulas <- c(
        lapply(seq(from = -1, to = 1, length = P), 
                function (rho) normalCopula(rho)),
        apply(merge(seq(from = -1, to = 1, length = P),
                        seq(from = 1, to = 30, length = P)), 1,
                function (p) tCopula(p[1], df = as.integer(p[2]), df.fixed = TRUE)),
        lapply(seq(from = .Machine$double.eps^0.5, to = 3, length = P),
                function (theta) claytonCopula(theta)),
        lapply(seq(from = 1, to = 5, length = P),
                function (theta) gumbelCopula(theta)),
        lapply(seq(from = -1, to = 1, length = P), 
                function (theta) fgmCopula(theta)),
        lapply(seq(from = 1, to = 5, length = P),
                function (theta) galambosCopula(theta)),       
        lapply(seq(from = -10, to = 10, length = P),
                function (theta) frankCopula(theta)),
        list(indepCopula()))

for (copula in copulas) {
    # Validate the h-function.
    u <- h(copula, XV[ , 1], XV[ , 2])
    uu <- vines:::hCopula(copula, XV[ , 1], XV[ , 2])
    stopifnot(isTRUE(all.equal(u, uu, T)))
    
    # Validate the inverse of the h-function.
    x <- hinverse(copula, u, XV[ , 2])
    xx <- vines:::hinverseCopula(copula, u, XV[ , 2])
    stopifnot(isTRUE(all.equal(x, xx, T)))
}


X <- seq(from = 0, to = 1, length = N)
V <- seq(from = 0, to = 1, length = N)
XV <- merge(X, V)

copulas <- c(
        lapply(seq(from = -1, to = 1, length = P), 
                function (rho) normalCopula(rho)),
        apply(merge(seq(from = -1, to = 1, length = P),
                        seq(from = 1, to = 30, length = P)), 1,
                function (p) tCopula(p[1], df = as.integer(p[2]), df.fixed = TRUE)),
        lapply(seq(from = .Machine$double.eps^0.5, to = 100, length = P),
                function (theta) claytonCopula(theta)),
        lapply(seq(from = 1, to = 100, length = P),
                function (theta) gumbelCopula(theta)),
        lapply(seq(from = -1, to = 1, length = P),
                function (theta) fgmCopula(theta)),
        lapply(seq(from = 1, to = 25, length = P),
                function (theta) galambosCopula(theta)),
        lapply(seq(from = -100, to = 100, length = P),
                function (theta) frankCopula(theta)),       
        list(indepCopula()))

for (copula in copulas) {
    # Check h(0, v) == 0.
    u <- h(copula, rep(0, N), V)
    stopifnot(isTRUE(all.equal(u, rep(0, N), T)))

    # Check h(1, v) == 1.
    u <- h(copula, rep(1, N), V)
    stopifnot(isTRUE(all.equal(u, rep(1, N), T)))

    # Check for finite return values in [0,1].
    u <- h(copula, XV[ , 1], XV[ , 2])
    stopifnot(all(is.finite(u)))
    stopifnot(all(u >= 0 & u <= 1))

    # Check for finite return values in [0,1].
    x <- hinverse(copula, u, XV[ , 2])
    stopifnot(all(is.finite(x)))
    stopifnot(all(x >= 0 & x <= 1))

    # Check hinverse(0, v) == 0.
    x <- hinverse(copula, rep(0, N), V)
    stopifnot(isTRUE(all.equal(x, rep(0, N), T)))

    # Check hinverse(1, v) == 1.
    x <- hinverse(copula, rep(1, N), V)
    stopifnot(isTRUE(all.equal(x, rep(1, N), T)))
}
