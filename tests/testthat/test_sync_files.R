test_that("sync_files_local and remove_local_files work", {
    expect_error(sync_files_local(), "have to be defined")
    tmpdir <- tempdir()
    dir.create(paste0(tmpdir, "/a"))
    dir.create(paste0(tmpdir, "/b"))
    fls_a <- file.create(paste0(tmpdir, "/a/", 1:5))
    expect_true(all(sync_files_local(1:5, paste0(tmpdir, "/a"),
                                     paste0(tmpdir, "/b"))))
    expect_true(length(sync_files_local(1:5, paste0(tmpdir, "/a"),
                                        paste0(tmpdir, "/b"))) == 0)
    expect_true(file.exists(paste0(tmpdir, "/b/3")))

    remove_local_files(1:4, paste0(tmpdir, "/b"))
    expect_false(file.exists(paste0(tmpdir, "/b/3")))
    expect_true(file.exists(paste0(tmpdir, "/b/5")))
})

