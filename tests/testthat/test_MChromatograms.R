test_that("c,MChromatograms works", {
    chr1 <- MSnbase::Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- MSnbase::Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- MSnbase::Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- MSnbase::MChromatograms(list(chr1, chr2, chr3))

    res <- .bind_rows_chromatograms(chrs)
    expect_equal(res, chrs)
    res <- .bind_rows_chromatograms(chrs, chrs)
    expect_true(nrow(res) == 6)
    expect_equal(res[3, ], chrs[3, ])
    expect_equal(res[6, ], chrs[3, ])

    res_2 <- c(chrs, chrs)
    expect_equal(res, res_2)

    chrs <- MSnbase::MChromatograms(list(chr1, chr2, chr3), ncol = 3)

    res <- c(chrs, chrs)
    expect_equal(chrs[1, 2], res[1, 2])
    expect_equal(chrs[1, 3], res[2, 3])

    chrs_2 <- MSnbase::MChromatograms(list(chr1, chr2), ncol = 2)
    expect_error(c(chrs, chrs_2), "must match")
})

test_that("plotOverlay works", {
    chr1 <- MSnbase::Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- MSnbase::Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- MSnbase::Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- MSnbase::MChromatograms(list(chr1, chr2, chr3))
    ## plotOverlay(chrs)
    ## plotOverlay(chrs, col = c("red", "blue", "green"))
})
