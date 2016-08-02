#' Cross-Validation Scheme Sampling
#'
#' \code{sample_cv} returns a matrix with cross-validation (CV) folds.
#'
#' @param x A vector with hybrid names.
#' @param n_mother A numeric or integer scalar specifying the size of the
#'  sampled maternal parents in the training set.
#' @param n_father A numeric or integer scalar specifying the size of the
#'  sampled paternal parents in the training set.
#' @param n_hyb_trn A numeric or integer scalar specifying the size of the
#'  training set.
#' @param min_size A numeric or integer scalar specifying the minimum size of
#'  each test set class (VS0, VS1, VS2).
#' @param rounds A numeric or integer speciyfing the number of cross-validation
#'  runs to sample.
#' @param hybrid_split A character scalar specifying the character that
#'  separates the maternal from the paternal parent in \code{x}.
#' @param progress A logical scalar. Set to 'TRUE' if a progress bar is desired.
#'
#' @return A matrix of dimensions \code{length(x)} by \code{rounds} with
#'  columns containing assignments of hybrids to either TRN, T0, T1 or T2.
#' @examples
#' data(hybrid_nms)
#' cv_mat <- sample_cv(hybrid_nms, n_mother = 39, n_father = 33,
#'                     n_hyb_trn = 200, min_size = 20, rounds = 20L,
#'                     hybrid_split = "_")
#' str(cv_mat)
#' @export
sample_cv <- function(x, n_father, n_mother, n_hyb_trn, min_size, rounds,
                      hybrid_split, progress = FALSE) {
  # Input testing.
  scalar_types <- vapply(list(n_mother, n_father, n_hyb_trn, min_size, rounds),
                         FUN = typeof, FUN.VALUE = character(1))
  stopifnot(all(scalar_types %in% c("double", "integer")))
  scalar_lengths <- vapply(list(n_mother, n_father, n_hyb_trn, min_size,
                                rounds),
                           FUN = length, FUN.VALUE = integer(1))
  stopifnot(all(scalar_lengths == 1))
  stopifnot(typeof(x) == "character")
  stopifnot(all(grepl(hybrid_split, x)))
  parent_number <- vapply(strsplit(x, split = hybrid_split), FUN = length,
                          FUN.VALUE = numeric(1))
  stopifnot(all(parent_number == 2))

  # Matrix for results
  container <- matrix(NA_character_, nrow = length(x), ncol = rounds)

  # Progress bar
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 1, max = rounds, initial = 1, char = "=",
                                style = 3, width = 80, label = min_size)
  }

  # Initial counter that will be incremented to a value of 'rounds' + 1.
  i <- 1L
  # Parents of hybrids in 'x'.
  mother <- vapply(strsplit(x, split = hybrid_split),
                   FUN = "[[", 1, FUN.VALUE = character(1))
  unique_mother <- unique(mother)
  father <- vapply(strsplit(x, split = hybrid_split),
                         FUN = "[[", 2, FUN.VALUE = character(1))
  unique_father <- unique(father)

  while (i != (rounds + 1)) {
    # Update the progress bar
    if (isTRUE(progress)) {
      utils::setTxtProgressBar(pb, value = i)
    }

    # Sample mother and father training set parents.
    m_trn <- sample(unique_mother, size = n_mother, replace = FALSE)
    f_trn <- sample(unique_father, size = n_father, replace = FALSE)
    # Assemble all hybrids whose parents are both tested.
    t2_trn <- x[mother %in% m_trn & father %in% f_trn]
    if (!all(c(m_trn, f_trn) %in% unlist(strsplit(t2_trn, split = "_")))) {
      next
    }
    t0 <- x[!mother %in% m_trn & !father %in% f_trn]
    if (any(unlist(strsplit(t0, split = "_")) %in% c(m_trn, f_trn))) {
      stop("Training set parents in T0 set")
    }
    n_t0 <- length(t0)
    if (n_t0 < min_size) {
      next
    }
    # Sample T2-hybrids.
    trn_mat <- matrix(paste(rep(m_trn, times = length(f_trn)),
                            rep(f_trn, each = length(m_trn)), sep = "_"),
                      nrow = length(m_trn), ncol = length(f_trn),
                      dimnames = list(m_trn, f_trn))
    trn_mat_copy <- trn_mat
    trn_mat[] <- ifelse(trn_mat %in% t2_trn, yes = "1", no = "0")
    t2_idx <- c()
    for (j in seq_len(prod(dim(trn_mat)))) {
      if (trn_mat[j] == "1") {
        trn_mat[j] <- "0"
        if (all(c(rowSums(trn_mat == "1"), colSums(trn_mat == "1")) != 0)) {
          t2_idx <- c(t2_idx, j)
        } else {
          trn_mat[j] <- "1"
        }
      }
      # Stop T2 sampling once the TRN is large enough.
      if ((length(t2_trn) - length(t2_idx)) == n_hyb_trn) break
    }
    t2 <- trn_mat_copy[t2_idx]
    n_t2 <- length(t2)
    if (n_t2 < min_size) {
      next
    }
    # Assign all non T2-hybrids to TRN.
    trn <- t2_trn[!t2_trn %in% t2]
    n_trn <- length(trn)
    if (n_trn < n_hyb_trn) {
      next
    }
    stopifnot(all(c(m_trn, f_trn) %in% unlist(strsplit(trn, split = "_"))))

    t1 <- setdiff(x, c(t0, t2, trn))
    # Check whether all T1-hybrids have exactly one parent that is part of the
    # training set.
    t1_mother <- vapply(strsplit(t1, split = "_"), "[[", 1,
                      FUN.VALUE = character(1))
    t1_father <- vapply(strsplit(t1, split = "_"), "[[", 2,
                       FUN.VALUE = character(1))
    t10 <- t1[t1_mother %in% m_trn & !t1_father %in% f_trn]
    t01 <- t1[!t1_mother %in% m_trn & t1_father %in% f_trn]
    stopifnot(identical(setdiff(t1, t10), t01))
    stopifnot(identical(setdiff(t1, t01), t10))
    n_t1 <- length(t1)
    if (n_t1 < min_size) {
      print("T1 too small")
      next
    }
    stopifnot(sum(n_t0, n_t1, n_t2, n_trn) == length(x))

    cv <- utils::stack(list(TRN = trn,
                            T0 = t0,
                            T1 = t1,
                            T2 = t2))
    cv <- cv[order(cv$values), ]
    container[, i] <- as.character(cv$ind)
    if (i == 1) {
      rownames(container) <- cv$values
    }

    # Check whether all TST-sets are unique for each run.
    if (i != 1) {
      # Check whether all TST-sets are unique for each partition.
      if (!all(vapply(c("T0", "T1", "T2", "TRN"), FUN = function(tst_set) {
        anyDuplicated(container[as.character(cv$ind) == tst_set, 1:i],
                      MARGIN = 2) == 0
      }, FUN.VALUE = logical(1)))) next
    }
    i <- i + 1L
  }
  rownames(container) <- paste0("DF_", rownames(container))
  out <- list(container)
  names(out) <-   paste0("ps",
                         100 * round(n_mother / length(unique_mother),
                                     digits = 2),
                         100 * round(n_father / length(unique_father),
                                     digits = 2),
                         "_trn=", n_hyb_trn, "_min_size=", min_size,
                         "_m=", n_mother, "_f=", n_father, ".RDS")
  out
}
