test_that("matchRtMz works", {
    x <- data.frame(mz = c(23.4, 45.6, 56.9, 76.5, 76.5, 76.5, 80.1),
                    rt = c(12, 34, 59, 34, 67, 65, 67))

    set.seed(123)
    y <- rbind(x, x)
    y$mz <- y$mz + rnorm(nrow(y), sd = 0.0002)
    y$rt[1:nrow(x)] <- x$rt + 2
    y <- y[order(y$mz), ]

    expect_error(matchRtMz(1:3, x), "have to be ")
    expect_error(matchRtMz(x[, 1, drop = FALSE], y), "Columns")
    expect_error(matchRtMz(x, y, mzcol = "other"), "Columns")

    res <- matchRtMz(x, y, ppm = 0)
    expect_true(all(is.na(res)))

    res <- matchRtMz(x, y)
    expect_true(length(res) == nrow(x))
    expect_true(is.integer(res))
    expect_equal(res, c(2L, 4L, 5L, 7L, 8L, 9L, 13L))

    res <- matchRtMz(x, y, duplicates = "keep")
    expect_true(is.list(res))
    expect_equal(res[[1]], c(1, 2))

    res <- matchRtMz(x, y, duplicates = "keep", ppm = 10)
    expect_equal(res[[1]], 2)

    res <- matchRtMz(x, y, ppm = 3)
    expect_true(is.na(res[[1]]))
})
