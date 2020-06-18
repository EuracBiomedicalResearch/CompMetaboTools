test_that("SimillarRtimeParam works", {
    res <- SimilarRtimeParam(4)
    expect_true(res@diffRt == 4)

    expect_error(SimilarRtimeParam(1:2), "positive numeric")
    expect_error(SimilarRtimeParam(-1), "positive numeric")

    expect_error(SimilarRtimeParam(method = "some"), "should be one")
    expect_error(new("SimilarRtimeParam", method = c("groupClosest", "greedy")),
                 "should be either")
    
    prm <- SimilarRtimeParam(3)
    
    expect_error(groupFeatures(xod, prm), "No feature definitions")
    res <- groupFeatures(xodg, prm)
    expect_true(any(colnames(featureDefinitions(res)) == "feature_group"))

    res <- groupFeatures(res, prm)

    tmp <- xodg
    featureDefinitions(tmp)$ms_level[c(1:3, 5)] <- 2
    res_2 <- groupFeatures(tmp, prm)
    fgs_2 <- featureDefinitions(res_2)$feature_group
    expect_true(all(is.na(fgs_2[c(1:3, 5)])))
    
    expect_error(groupFeatures(res, prm, msLevel = 1:2), "Currently only")
})

test_that("EicCorrelationParam works", {
    res <- EicCorrelationParam(threshold = 4)
    expect_equal(res@threshold, 4)

    expect_error(EicCorrelationParam(threshold = c(1, 2)), "positive numeric")
    expect_error(EicCorrelationParam(n = 1:2), "positive numeric")
    expect_error(EicCorrelationParam(clean = c(TRUE, FALSE)), "length 1")
    expect_error(EicCorrelationParam(value = "other"))

    ## n bigger than 3
    expect_error(groupFeatures(xodg, param = EicCorrelationParam(n = 5)),
                 "smaller or")
    ## no feature definitions
    expect_error(groupFeatures(xod, param = EicCorrelationParam()), "No")
    ## MS level length > 1
    expect_error(
        groupFeatures(xodg, param = EicCorrelationParam(), msLevel = 1:2),
        "single MS level")

    idx <- c(9, 30, 32, 45, 61, 99, 100, 104, 110, 115, 120, 121, 122, 123)
    tmp <- xodg
    featureDefinitions(tmp)$feature_group <- NA
    featureDefinitions(tmp)$feature_group[idx] <- "FG"
    res <- groupFeatures(tmp, param = EicCorrelationParam())
    expect_true(all(is.na(featureDefinitions(res)$feature_group[-idx])))

    expect_true(length(unique(featureDefinitions(res)$feature_group)) < length(idx))

    featureDefinitions(tmp)$feature_group <- NULL
    featureDefinitions(tmp)$ms_level[idx] <- 2

    res_2 <- groupFeatures(tmp, param = EicCorrelationParam(), msLevel = 2)
    expect_equal(featureDefinitions(res)$feature_group,
                 featureDefinitions(res_2)$feature_group)
})

test_that("AbundanceCorrelationParam works", {
    prm <- AbundanceCorrelationParam(threshold = 0.5, value = "maxo")
    expect_equal(prm@threshold, 0.5)
    expect_equal(prm@value, "maxo")
    expect_equal(prm@method, "maxint")

    expect_error(AbundanceCorrelationParam(subset = "4"), "integer")
    
    expect_error(groupFeatures(xod, AbundanceCorrelationParam()), "feature")
    expect_error(
        groupFeatures(xodg, AbundanceCorrelationParam(subset = c(1, 4, 5))),
        "has to be between")
    expect_error(
        groupFeatures(xodg, AbundanceCorrelationParam(subset = c(TRUE))),
        "Can not calculate")

    res <- groupFeatures(xodg, AbundanceCorrelationParam())
    expect_true(any(colnames(featureDefinitions(res)) == "feature_group"))
    expect_true(length(unique(featureDefinitions(res)$feature_group)) <
                nrow(featureDefinitions(res)))
    res_2 <- groupFeatures(xodg, AbundanceCorrelationParam(subset = c(2, 3)))

    plotFeatureGroups(res_2)
    expect_error(plotFeatureGroups(res_2, featureGroups = "a"), "None of the")
    expect_error(plotFeatureGroups(xodg), "No feature groups")
    
    ## With pre-defined grps.
    tmp <- xodg
    featureDefinitions(tmp)$feature_group <- "FG.2"
    idx <- c(4, 12, 23, 56)
    featureDefinitions(tmp)$ms_level[idx] <- 2

    res <- groupFeatures(tmp, AbundanceCorrelationParam(), msLevel = 1)
    expect_true(all(featureGroups(res)[idx] == "FG.2"))
    expect_true(all(featureGroups(res)[-idx] != "FG.2"))
    res_2 <- groupFeatures(tmp, AbundanceCorrelationParam(), msLevel = 2)
    expect_true(all(featureGroups(res_2)[-idx] == "FG.2"))
    expect_true(all(featureGroups(res_2)[idx] != "FG.2"))

})

test_that(".format_groups works", {
    res <- .format_groups(1:3)
    expect_equal(res, c("1", "2", "3"))
    res <- .format_groups(1:10)
    expect_equal(res, c("01", "02", "03", "04", "05", "06", "07", "08",
                        "09", "10"))
})

test_that("featureGroupPseudoSpectrum works", {
    fvals <- featureValues(xodgg, method = "maxint", value = "maxo")
    ## 3 feature
    ft_idx <- which(featureGroups(xodgg) == "FG.009.1")
    res <- featureGroupPseudoSpectrum("FG.009.1", xodgg, fvals = fvals,
                                      intensityFun = median)
    expect_true(is(res, "Spectrum1"))
    expect_true(peaksCount(res) == 3)
    expect_true(validObject(res))
    expect_equal(intensity(res), apply(fvals[ft_idx, ], MARGIN = 1,
                                       median, na.rm = TRUE))
    expect_equal(mz(res), featureDefinitions(xodgg)$mzmed[ft_idx])
    expect_equal(rtime(res), median(featureDefinitions(xodgg)$rtmed[ft_idx]))
    
    ## 1 feature
    res <- featureGroupPseudoSpectrum("FG.009.2", xodgg, fvals = fvals,
                                      intensityFun = median)
    ft_idx <- which(featureGroups(xodgg) == "FG.009.2")
    expect_true(is(res, "Spectrum1"))
    expect_true(peaksCount(res) == 1)
    expect_true(validObject(res))
    expect_equal(unname(intensity(res)),
                 unname(median(fvals[ft_idx, ], na.rm = TRUE)))
    expect_equal(mz(res), featureDefinitions(xodgg)$mzmed[ft_idx])
    expect_equal(rtime(res), median(featureDefinitions(xodgg)$rtmed[ft_idx]))

    expect_error(
        featureGroupPseudoSpectrum("FG.009.1", xodgg, fvals = fvals, n = 12),
        "has to be an integer")
})

test_that("featureGroupFullScan works", {
    fvals <- featureValues(xodgg, method = "maxint", value = "maxo")
    ## 3 feature
    res <- featureGroupFullScan("FG.009.1", xodgg, fvals = fvals)
    ft_idx <- which(featureGroups(xodgg) == "FG.009.1")
    expect_true(is(res, "Spectrum1"))
    expect_true(
        abs(rtime(res) -
            median(featureDefinitions(xodgg)[ft_idx, "rtmed"])) < 0.1)
    expect_true(all(featureDefinitions(xodgg)[ft_idx, "mzmed"] %in% mz(res)))

    ## 1 feature
    res <- featureGroupFullScan("FG.009.2", xodgg, fvals = fvals)
    ft_idx <- which(featureGroups(xodgg) == "FG.009.2")
    expect_true(is(res, "Spectrum1"))
    expect_true(
        abs(rtime(res) -
            median(featureDefinitions(xodgg)[ft_idx, "rtmed"])) < 0.8)
    expect_true(all(featureDefinitions(xodgg)[ft_idx, "mzmed"] %in% mz(res)))
})

test_that("featureGroupSpectra works", {
    ## Errors
    expect_error(featureGroupSpectra(xodgg, subset = 1:5), "an integer")
    expect_error(featureGroupSpectra(xod), "feature definitions")
    expect_error(featureGroupSpectra(xodgg, featureGroup = c("a")), "all feature")
    
    ## Get all of them
    res_all <- featureGroupSpectra(xodgg)
    expect_true(is(res_all, "Spectra"))
    expect_equal(mcols(res_all)$feature_group, unique(featureGroups(xodgg)))
    expect_equal(unname(peaksCount(res_all)),
                 unname(lengths(mcols(res_all)$feature_id)))
    
    ## Get them in a subset
    res_sub <- featureGroupSpectra(xodgg, subset = c(1, 3))
    expect_true(sum(is.na(rtime(res_sub))) == 59)
    
    ## Get only selected ones
    res <- featureGroupSpectra(xodgg, featureGroup = c("FG.009.1", "FG.009.2"))
    expect_true(length(res) == 2)
    expect_equal(mcols(res)$feature_group, c("FG.009.1", "FG.009.2"))
    idx <- which(mcols(res_all)$feature_group %in% c("FG.009.1", "FG.009.2"))
    expect_equal(res[[1]], res_all[[idx[1]]])
    expect_equal(res[[2]], res_all[[idx[2]]])
})
