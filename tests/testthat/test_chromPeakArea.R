test_that("chromPeakArea works", {
    rts <- c(2500, 2700)
    mzs <- cbind(c(400, 300), c(500, 400))
    res <- chromPeakArea(xod, rt = rts, diffRt = 100, mz = mzs)
    pks <- chromPeaks(xod)

    keep <- pks[, "mz"] > 400 & pks[, "mz"] < 500 &
        pks[, "rt"] > 2400 & pks[, "rt"] < 2600
    mzr <- range(pks[keep, c("mzmin", "mzmax")])
    expect_equal(mzr, unname(res[1, 1:2]))
    rtr <- range(pks[keep, c("rtmin", "rtmax")])
    expect_equal(rtr, unname(res[1, 3:4]))

    keep <- pks[, "mz"] > 300 & pks[, "mz"] < 400 &
        pks[, "rt"] > 2600 & pks[, "rt"] < 2800
    mzr <- range(pks[keep, c("mzmin", "mzmax")])
    expect_equal(mzr, unname(res[2, 1:2]))
    rtr <- range(pks[keep, c("rtmin", "rtmax")])
    expect_equal(rtr, unname(res[2, 3:4]))
    
    expect_error(chromPeakArea(xod, peakId = "some"), "implemented")
    expect_error(chromPeakArea(xod, rt = 2), "have to match")
})
