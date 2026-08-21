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
