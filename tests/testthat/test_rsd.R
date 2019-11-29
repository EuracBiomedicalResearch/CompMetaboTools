test_that("rsd, rowRsd works", {
    a <- c(5.3, 5.7, 2.5, 2.5)
    expect_equal(rsd(a), sd(a) / mean(a))

    res <- rowRsd(rbind(a, a))
    expect_equal(unname(res[1]), rsd(a))
})
