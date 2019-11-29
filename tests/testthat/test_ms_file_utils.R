test_that("extract_time_stamp works", {
    fl <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML",
                      package = "msdata")
    res <- extract_time_stamp(fl)
    expect_true(is(res, "POSIXt"))
})
