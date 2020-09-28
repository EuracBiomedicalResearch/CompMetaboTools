test_that("joyPlot works", {
    library(xcms)
    data(faahko_sub)
    dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
    one <- filterFile(faahko_sub, 1L)

    joyPlot(one, rt = c(4000, 4200))
    
    mzs <- seq(400, 430, by = 1)
    chr <- chromatogram(one, mz = cbind(mzs[-length(mzs)], mzs[-1]),
                        rt = c(4000, 4200), include = "none")
    joyPlot(chr)
    .joy_plot_chromatograms(chr, col = heat.colors(nrow(chr)))
})
