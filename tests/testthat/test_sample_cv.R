data("hybrid_nms")
# -----------------------------------------------------------------------------
context("Uniqueness of CV runs")
test_that("no CV runs is duplicated", {
  expect_identical(check_cv(sample_cv(hybrid_nms, n_mother = 40,
                                      n_father = 33,
                                      n_hyb_trn = 200,
                                      min_size = 25, rounds = 100,
                                      hybrid_split = "_")),
                   "success")
})
# -----------------------------------------------------------------------------
# EOF
