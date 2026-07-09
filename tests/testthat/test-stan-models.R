test_that(".stan_source_path() finds all shipped .stan files", {
  for (nm in c("model_betadir", "model_pureNMF", "model_test")) {
    path <- ClassTopics:::.stan_source_path(nm)
    expect_true(file.exists(path))
    expect_true(grepl(paste0(nm, "\\.stan$"), path))
  }
})

test_that(".stan_source_path() rejects unknown model names", {
  expect_error(ClassTopics:::.stan_source_path("not_a_real_model"))
})

test_that(".stan_cache_dir() returns a writable, existing directory", {
  cache_dir <- ClassTopics:::.stan_cache_dir()
  expect_true(dir.exists(cache_dir))
  # Confirm it is actually writable
  probe <- file.path(cache_dir, "ClassTopics_write_probe.tmp")
  on.exit(unlink(probe), add = TRUE)
  expect_true(file.create(probe))
})
