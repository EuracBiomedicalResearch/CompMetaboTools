library("testthat")
library("CompMetaboTools")
library(MsFeatures)

library(faahKO)

fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
         system.file('cdf/KO/ko16.CDF', package = "faahKO"),
         system.file('cdf/WT/wt19.CDF', package = "faahKO"))

od <- readMSData(fls, mode = "onDisk")
xod <- findChromPeaks(
    od, param = CentWaveParam(noise = 10000, snthresh = 40,
                              prefilter = c(3, 10000)))
pdp <- PeakDensityParam(sampleGroups = c(1, 1, 2))
xodg <- groupChromPeaks(xod, param = pdp)
xodgg <- groupFeatures(xodg, param = SimilarRtimeParam(4))
xodgg <- groupFeatures(xodgg, param = AbundanceSimilarityParam(threshold = 0.3))

test_check("CompMetaboTools")
