data("hybrid_nms")
# -----------------------------------------------------------------------------
context("Timer functionality")
test_that("job stops after 'max_time' minutes", {
expect_match(sample_cv(hybrid_nms, n_mother = 39, n_father = 33,
                       n_hybrid = 200,
                       min_size = 20, rounds = 1000L, hybrid_split = "_",
                       max_time = 0.1),
             regexp = "job took too long")
})

context("Uniqueness of CV runs")
test_that("no CV runs is duplicated", {
  expect_false(any(duplicated(sample_cv(hybrid_nms, n_mother = 39,
                                        n_father = 33,
                                        n_hybrid = 200,
                                        min_size = 20, rounds = 20L,
                                        hybrid_split = "_"), MARGIN = 2)))
})

# -----------------------------------------------------------------------------
context("Number of parents in TRN")
test_that("actual number of TRN-parents equal to specified number of parents",{
  expect_equal(sspredr::check_trn(sample_cv(hybrid_nms,
                                            n_mother = 39,
                                            n_father = 33,
                                            n_hybrid = 200,
                                            min_size = 20,
                                            rounds = 20L,
                                            hybrid_split = "_")),
               c(n_mother = 39, n_father = 33))
})


# -----------------------------------------------------------------------------
# EOF
