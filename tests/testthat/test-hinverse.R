context("hinverse")

test_that("hinverse() returns finite values in [0,1]", {
    for (copula in copulas) {
        x <- hinverse(copula, XV[ , 1], XV[ , 2])
        expect_true(all(is.finite(x)))
        expect_true(all(x >= 0 & x <= 1))
    }
})

test_that("hinverse() computes the correct values", {
    XV <- subset(XV, X >= 0.25 & X <= 0.75 & V >= 0.25 & V <= 0.75)
    for (copula in copulas) {
        x <- hinverse(copula, XV[ , 1], XV[ , 2])
        xx <- vines:::hinverseCopula(copula, XV[ , 1], XV[ , 2],
                                     eps = .Machine$double.eps^0.5)
        expect_equal(x, xx, tolerance = tol)
    }
})

test_that("hinverse(0, v) == 0", {
    for (copula in copulas) {
        x <- hinverse(copula, rep(0, n), V)
        expect_equal(x, rep(0, n), tolerance = tol)
    }
})

test_that("hinverse(1, v) == 1", {
    for (copula in copulas) {
        x <- hinverse(copula, rep(1, n), V)
        expect_equal(x, rep(1, n), tolerance = tol)
    }
})
