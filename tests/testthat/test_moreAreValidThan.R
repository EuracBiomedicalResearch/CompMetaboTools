test_that("moreAreValidThan works", {
    x <- rbind(
        c(NA, 3, 4, 1, 3, NA, 4, NA),
        c(4, 2, 3, 4, 5, 5, 2, NA),
        c(NA, NA, NA, NA, NA, 3, 4, 5))

    res <- moreAreValidThan(x, prop = 0.5)
    expect_equal(res, c(TRUE, TRUE, FALSE))

    res <- moreAreValidThan(x, prop = 0.5, f = c(1, 1, 1, 1, 1, 2, 2, 2))
    expect_equal(res, c(TRUE, TRUE, TRUE))

    res <- moreAreValidThan(x, prop = 0.5, f = c(1, 1, 1, 1, 1, 2, 2, 2),
                            condition = all)
    expect_equal(res, c(FALSE, TRUE, FALSE))

    expect_error(moreAreValidThan(x, f = 1:3))
})
