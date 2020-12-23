test_that(".group_logic_matrix works", {
    xmat <- rbind(c(TRUE, FALSE, FALSE, FALSE),
                  c(FALSE, TRUE, FALSE, FALSE),
                  c(FALSE, FALSE, TRUE, FALSE),
                  c(FALSE, FALSE, FALSE, TRUE))
    expect_error(.group_logic_matrix(xmat[1:3, ]))
    res <- .group_logic_matrix(xmat)
    expect_true(length(res) == nrow(xmat))
    expect_equal(res, list(1, 2, 3, 4))

    xmat <- rbind(c(TRUE, FALSE, FALSE, FALSE, TRUE),
                  c(FALSE, TRUE, FALSE, FALSE, FALSE),
                  c(FALSE, FALSE, TRUE, FALSE, TRUE),
                  c(FALSE, FALSE, FALSE, TRUE, FALSE),
                  c(TRUE, FALSE, TRUE, FALSE, TRUE))
    res <- CompMetaboTools:::.group_logic_matrix(xmat)
    expect_equal(res, list(c(1, 3, 5), 2, 4))

    xcor <- matrix(FALSE, ncol = 13, nrow = 13)
    for (i in 1:13)
        xcor[i, i] <- TRUE
    xcor[8, 6] <- TRUE
    xcor[8, 7] <- TRUE
    xcor[9, 7] <- TRUE
    xcor[11, 7] <- TRUE
    xcor[6, 8] <- TRUE
    xcor[7, 8] <- TRUE
    xcor[10, 8] <- TRUE
    xcor[13, 8] <- TRUE
    xcor[7, 9] <- TRUE
    xcor[8, 10] <- TRUE
    xcor[7, 11] <- TRUE
    xcor[12, 11] <- TRUE
    xcor[11, 12] <- TRUE
    xcor[8, 13] <- TRUE
    res <- .group_logic_matrix(xcor)
    expect_equal(res, list(1, 2, 3, 4, 5, c(6:13)))
    
    xcor <- matrix(FALSE, ncol = 10, nrow = 10)
    for (i in seq_len(ncol(xcor))) {
        xcor[i, i] <- TRUE
    }
    xcor[1, 4] <- TRUE
    xcor[4, 1] <- TRUE
    xcor[2, 8] <- TRUE
    xcor[8, 2] <- TRUE
    xcor[3, 9] <- TRUE
    xcor[9, 3] <- TRUE
    xcor[8, 9] <- TRUE
    xcor[9, 8] <- TRUE
    res <- .group_logic_matrix(xcor)
    expect_equal(res, list(c(1, 4), c(2, 3, 8, 9), 5, 6, 7, 10))
})

test_that(".index_list_to_factor works", {
    x <- list(c(1, 5, 2), c(3, 4), c(6), 7)
    res <- .index_list_to_factor(x)
    expect_equal(res, factor(c(1, 1, 2, 2, 1, 3, 4)))
})

test_that("groupByCorrelation works", {
    x <- rbind(c(1, 2, 3, 4),
               c(2, 4, 6, 8),
               c(0, 2, 1, 2),
               c(1, 3, 4, 5))
    res <- groupByCorrelation(x)
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1, 1, 2, 1)))

    res_2 <- groupByCorrelation(x, inclusive = TRUE)
    expect_equal(res, res_2)
    
    expect_error(groupByCorrelation(x, threshold = c(0.4, 0.3)), "length 1")

    x <- rbind(x,
               c(2, 4, 6, 9),
               c(1, 4, 1, 4),
               c(1, 2, 3, 4))
    f <- c(1, 2, 2, 1, 1, 2, 2)
    res <- groupByCorrelation(x, f = f)
    expect_equal(res, factor(c("1.001", "2.001", "2.002", "1.001",
                               "1.001", "2.002", "2.001")))

    f <- c(1, 2, NA, NA, 1, 2, 2)
    res <- groupByCorrelation(x, f = f)
    expect_equal(res, factor(c("1.001", "2.001", NA, NA, "1.001",
                               "2.002", "2.001")))
    
    expect_error(groupByCorrelation(x, f = 3), "its length has to ")
})

test_that("groupEicCorrelation works", {
    set.seed(123)
    chr1 <- MSnbase::Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- MSnbase::Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- MSnbase::Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- MSnbase::MChromatograms(list(chr1, chr2, chr3))

    res <- groupEicCorrelation(chrs, align = "closest")
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1L, 2L, 1L)))
    res <- groupEicCorrelation(chrs, align = "closest", tolerance = 0)
    expect_equal(res, factor(c(1L, 2L, 3L)))
    
    chrs <- MSnbase::MChromatograms(list(chr1, chr2, chr3, chr1, chr2, chr3,
                                        chr2, chr3, chr1), ncol = 3)
    res <- groupEicCorrelation(chrs)
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1L, 2L, 3L)))

    res <- groupEicCorrelation(chrs, aggregationFun = max, inclusive = TRUE,
                               align = "closest")
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1L, 1L, 1L)))

    res <- groupEicCorrelation(chrs, aggregationFun = max, inclusive = FALSE,
                               align = "closest")
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1L, 2L, 1L)))

    res <- groupEicCorrelation(chrs, aggregationFun = median,
                               align = "closest")
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1L, 2L, 1L)))
})

test_that("groupToSinglePolarityPairs works", {
    x <- rbind(
        c(4, 3, 5, 1),
        c(4, 2, 5, 1),
        c(4, 3, 4, 1),
        c(4, 3, 4, 1),
        c(4, 4, 4, 9),
        c(4, 4, 4, 9),
        c(4, 4, 4, 9))

    expect_error(groupToSinglePolarityPairs(f = 1:4), "are required")
    expect_error(groupToSinglePolarityPairs(f = rep(1, nrow(x)), polarity = 4,
                                fvals = x), "has to match")
    expect_error(groupToSinglePolarityPairs(f = rep(1, nrow(x)),
                                            polarity = rep(1, nrow(x)),
                                            fvals = x[1:4, ]), "has to match")

    res <- groupToSinglePolarityPairs(f = rep(1, nrow(x)),
                                      polarity = rep(1, nrow(x)),
                                      x)
    expect_true(length(res) == length(levels(res)))

    pol <- c("POS", "NEG", "POS", "NEG", "POS", "NEG", "POS")
    res <- groupToSinglePolarityPairs(f = c(1, 1, 1, 1, 1, 2, 2),
                                      polarity = pol, x)
    expect_equal(res, factor(c("1.2", "1.2", "1.1", "1.1", "1.3", "2.1", "2.1")))
})
