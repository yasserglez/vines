n <- 7  # Number of values of each variable
np <- 5  # Number of values of each parameter
tol <- 0.01  # Ignore differences smaller than tol

X <- seq(from = 0, to = 1, length = n)
V <- seq(from = 0, to = 1, length = n)
XV <- merge(X, V)
colnames(XV) <- c("X", "V")

copulas <- c(
    lapply(seq(from = -1, to = 1, length = np),
           function (rho) normalCopula(rho)),
    apply(merge(seq(from = -1, to = 1, length = np),
                seq(from = 1, to = 30, length = np)), 1,
          function (p) tCopula(p[1], df = as.integer(p[2]), df.fixed = TRUE)),
    lapply(seq(from = .Machine$double.eps^0.5, to = 50, length = np),
           function (theta) claytonCopula(theta, use.indepC = "FALSE")),
    lapply(seq(from = 1, to = 50, length = np),
           function (theta) gumbelCopula(theta, use.indepC = "FALSE")),
    lapply(seq(from = -1, to = 1, length = np),
           function (theta) fgmCopula(theta)),
    lapply(seq(from = 1, to = 25, length = np),
           function (theta) galambosCopula(theta)),
    lapply(seq(from = -50, to = 50, length = np),
           function (theta) frankCopula(theta, use.indepC = "FALSE")),
    list(indepCopula()))
