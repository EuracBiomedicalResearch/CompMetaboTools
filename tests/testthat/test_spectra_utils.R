test_that("extractSpectraData works", {
    fl <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML",
                      package = "msdata")
    data <- MSnbase::filterRt(MSnbase::readMSData(fl, mode = "onDisk"),
                              rt = c(1, 6))
    sps <- MSnbase::spectra(data)

    res <- extractSpectraData(sps)
    expect_true(is(res, "DataFrame"))
    expect_true(all(c("mz", "intensity") %in% colnames(res)))
    expect_true(is(res$mz, "NumericList"))
    expect_true(is(res$intensity, "NumericList"))

    res <- Spectra::Spectra(res)
    expect_true(is(res, "Spectra"))
    expect_equal(Spectra::msLevel(res), unname(MSnbase::msLevel(data)))

    expect_error(extractSpectraData(1:10), "should be either a 'list'")
    
    res <- extractSpectraData(data)
    expect_true(is(res, "DataFrame"))
    expect_true(all(c("mz", "intensity") %in% colnames(res)))
    expect_true(is(res$mz, "NumericList"))
    expect_true(is(res$intensity, "NumericList"))

    res <- Spectra::Spectra(res)
    expect_true(is(res, "Spectra"))
    expect_equal(Spectra::msLevel(res), unname(MSnbase::msLevel(data)))

    spctra <- MSnbase::MSpectra(sps)
    mcols(spctra)$new_col <- "a"
    res <- extractSpectraData(sps)
    expect_true(all(res$new_col == "a"))
})
