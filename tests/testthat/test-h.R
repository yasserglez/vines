context("h")

test_that("h() returns finite values in [0,1]", {
    for (copula in copulas) {
        u <- h(copula, XV[ , 1], XV[ , 2])
        expect_true(all(is.finite(u)))
        expect_true(all(u >= 0 & u <= 1))
    }
})

test_that("h() computes the correct values", {
    XV <- subset(XV, X >= 0.25 & X <= 0.75 & V >= 0.25 & V <= 0.75)
    for (copula in copulas) {
        u <- h(copula, XV[ , 1], XV[ , 2])
        uu <- vines:::hCopula(copula, XV[ , 1], XV[ , 2],
                              eps = .Machine$double.eps^0.5)
        expect_equal(u, uu, tolerance = tol)
    }
})

test_that("h(0, v) == 0", {
    for (copula in copulas) {
        u <- h(copula, rep(0, n), V)
        expect_equal(u, rep(0, n), tolerance = tol)
    }
})

test_that("h(1, v) == 1", {
    for (copula in copulas) {
        u <- h(copula, rep(1, n), V)
        expect_equal(u, rep(1, n), tolerance = tol)
    }
})
