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

    expect_warning(res <- groupFeatures(res, prm), "Found existing")

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
    featureDefinitions(tmp)$feature_group[idx] <- "FG.1"
    res <- groupFeatures(tmp, param = EicCorrelationParam())
    expect_true(all(is.na(featureDefinitions(res)$feature_group[-idx])))

    expect_true(length(unique(featureDefinitions(res)$feature_group)) < length(idx))

    featureDefinitions(tmp)$feature_group <- NULL
    featureDefinitions(tmp)$ms_level[idx] <- 2

    res_2 <- groupFeatures(tmp, param = EicCorrelationParam(), msLevel = 2)
    expect_equal(featureDefinitions(res)$feature_group,
                 featureDefinitions(res_2)$feature_group)
})
