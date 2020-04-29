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
    res <- .group_logic_matrix(xmat)
    expect_equal(res, list(c(1, 5, 3), 2, 4))
})

test_that("groupByCorrelation works", {
    x <- rbind(c(1, 2, 3, 4),
               c(2, 4, 6, 8),
               c(0, 2, 1, 2),
               c(1, 3, 4, 5))
    res <- groupByCorrelation(x)
    expect_equal(res, list(c(1, 2, 4), 3))

    expect_error(groupByCorrelation(x, threshold = c(0.4, 0.3)), "length 1")
})
