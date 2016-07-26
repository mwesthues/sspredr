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
test_cv <- function(x, n_father, n_mother, hybrid_split = "_") {
  dat <- as.data.frame(x)
  dat_cols <- ncol(dat)
  colnames(dat) <- paste0("Run_", seq_len(ncol(x)))
  dat$Hybrid <- rownames(x)
  rownames(dat) <- NULL

  split_id <- function(id, element, hybrid_split) {
  vapply(strsplit(id, split = hybrid_split), FUN = "[[", element,
         FUN.VALUE = character(1))
  }

  run_lst <- lapply(seq_len(dat_cols), FUN = function(i) {
    pick_col <- paste0("Run_", i)
    dat <- dat[, c(pick_col, "Hybrid")]
    ts_hybrids <- dat[dat[[get("pick_col")]] == "TS", "Hybrid"]
    split_lst <- list(n_mother = split_id(ts_hybrids, element = 1,
                                          hybrid_split = hybrid_split),
                      n_father = split_id(ts_hybrids, element = 2,
                                          hybrid_split = hybrid_split))
    ts_grp_length <- vapply(split_lst, FUN = function(x) length(unique(x)),
                            FUN.VALUE = integer(1))
  })
  run_mat <- do.call(rbind, run_lst)
  actual_ts <- vapply(seq_len(ncol(run_mat)),
                      FUN = function(i) unique(run_mat[, i]),
                      FUN.VALUE = integer(1))
  names(actual_ts) <- colnames(run_mat)
  actual_ts
}
test_that("actual number of TRN-parents equal to specified number of parents",{
  expect_equal(test_cv(sample_cv(hybrid_nms, n_mother = 39, n_father = 33,
                                 n_hybrid = 200, min_size = 20,
                                 rounds = 20L, hybrid_split = "_")),
               c(n_mother = 39, n_father = 33))
})


# -----------------------------------------------------------------------------
# EOF
