#' Check TRN Assembly
#' @keywords internal
#' @param x A matrix containing cross-validation folds inside columns and
#'  hybrid names in rows.
#' @param hybrid_split Character that separates parents in hybrid names.
#' @return Integer with the number of maternal and paternal training set lines.
check_trn <- function(x, hybrid_split = "_") {
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
