#' Cross-Validation Scheme Sampling
#'
#' \code{sample_cv} returns a matrix with cross-validation (CV) folds.
#'
#' @param x A vector with hybrid names.
#' @param n_mother A numeric or integer scalar specifying the size of the
#'  sampled maternal parents in the training set.
#' @param n_father A numeric or integer scalar specifying the size of the
#'  sampled paternal parents in the training set.
#' @param n_hybrid A numeric or integer scalar specifying the size of the
#'  training set.
#' @param min_size A numeric or integer scalar specifying the minimum size of
#'  each test set class (VS0, VS1, VS2).
#' @param rounds A numeric or integer speciyfing the number of cross-validation
#'  runs to sample.
#' @param hybrid_split A character scalar specifying the character that
#'  separates the maternal from the paternal parent in \code{x}.
#' @param max_time A numeric or integer scalar specifying the number of minutes
#'  after which the sampling procedure should be terminated.
#' @param progress A logical scalar. Set to 'TRUE' if a progress bar is desired.
#'
#' @return A matrix of dimensions \code{length(x)} by \code{rounds} with
#'  columns containing assignments of hybrids to either TS, VS0, VS1 or VS2.
#' @examples
#' data(hybrid_nms)
#' cv_mat <- sample_cv(hybrid_nms, n_mother = 39, n_father = 33, n_hybrid = 200,
#'                     min_size = 20, rounds = 20L, hybrid_split = "_")
#' dim(cv_mat)
#' cv_mat[1:5, 1:5]
#' @export
sample_cv <- function(x, n_father, n_mother, n_hybrid, min_size, rounds,
                      hybrid_split, max_time = NULL, progress = NULL) {

  # Input testing.
  scalar_types <- vapply(list(n_mother, n_father, n_hybrid, min_size, rounds),
                         FUN = typeof, FUN.VALUE = character(1))
  stopifnot(all(scalar_types %in% c("double", "integer")))
  if (!is.null(max_time)) {
    stopifnot(typeof(max_time) %in% c("double", "integer"))
    stopifnot(length(max_time) == 1)
  }
  scalar_lengths <- vapply(list(n_mother, n_father, n_hybrid, min_size, rounds),
                           FUN = length, FUN.VALUE = integer(1))
  stopifnot(all(scalar_lengths == 1))
  stopifnot(typeof(x) == "character")
  stopifnot(all(grepl(hybrid_split, x)))
  parent_number <- vapply(strsplit(x, split = hybrid_split), FUN = length,
                          FUN.VALUE = numeric(1))
  stopifnot(all(parent_number == 2))


  # Get the current system time as a reference. Some CV-assemblies may take an
  # excessive amount of time in bordeline cases that should be best avoided.
  # Hence, this timer will be used to terminate such jobs after 'max_time'
  # minutes.
  if (!is.null(max_time)) {
    init_time <- Sys.time()
  }
  # Matrix for results
  container <- matrix(NA_character_, nrow = length(x), ncol = rounds)
  # Progress bar
  if (!is.null(progress)) {
    pb <- utils::txtProgressBar(min = 1, max = rounds, initial = 1, char = "=",
                                style = 3, width = 80, label = min_size)
  }
  # Initial counter that will be incremented to a value of 'rounds' + 1.
  i <- 1
  while (i < (rounds + 1)) {
    # Terminate the job if it took longer than tolerated.
    if (!is.null(max_time)) {
      time_diff <- as.numeric(Sys.time() - init_time, units = "mins")
      if (time_diff >= max_time) {
        container <- "job took too long"
        break
      }
    }
    if (!is.null(progress)) {
      utils::setTxtProgressBar(pb, value = i)
    }

    # Sample one cross-validation run
    # The success of the sampling of the training set is dependent upon the
    # composition of the random subsets built from Dent (D_T) and Flint (F_T),
    # respectively.
    # Therefore, two while loops are run that break if the composition of either
    # of D_T or D_F is inadequate.
    dent <- unique(vapply(strsplit(x, split = hybrid_split),
                          FUN = "[[", 1, FUN.VALUE = character(1)))
    flint <- unique(vapply(strsplit(x, split = hybrid_split),
                           FUN = "[[", 2, FUN.VALUE = character(1)))
    n_t0 <- n_t1 <- n_t2 <- n_trn <- 0
    while_cond <- "failed"
    while (while_cond == "failed" ||
           any(c(n_t0, n_t1, n_t2, n_trn) < min_size)) {
      # Generate the random subsets of Dent and Flint genotypes
      dent_trn <- sample(x = dent, size = n_mother, replace = FALSE)
      flint_trn <- sample(x = flint, size = n_father, replace = FALSE)
      df_trn <- c(dent_trn, flint_trn)
      hybrid_list <- strsplit(x = x, split = hybrid_split)

      # Sample a random subset of n_hybrid training set hybrids from all hybrids
      # for which both the Dent and Flint parents were elements of dt and ft,
      # respectively.
      t2_index <- lapply(X = seq_len(length(hybrid_list)), FUN = function(i) {
          ind <- all(hybrid_list[[i]] %in% c(dent_trn, flint_trn))
          })
      t2_index <- unlist(t2_index)
      t2_hybrids <- x[t2_index]
      trn <- sample(x = t2_hybrids, size = n_hybrid, replace = FALSE)

      timer <- proc.time()["elapsed"]
      while_cond <- "worked"
      # Another constraint is necessary, namely that all lines in D_T and F_T
      # were parents of at least one hybrid in the training set.
      while (!all(grepl(paste0(unlist(strsplit(trn, split = "_")),
                               collapse = "|"), x = df_trn))) {
        trn <- sample(x = t2_hybrids, size = n_hybrid, replace = FALSE)
        # Break out of the loop if the given condition cannot be met due to
        # the composition of D_T or F_T, respectively.
          if ((proc.time()["elapsed"] - timer) > 1) {
            while_cond <- "failed"
            break
          }
        }

      # Hybrids, for which obth the Dent and the Flint parents were elements of
      # dt and ft, but were not elements of t2_train were assigned to the T2
      # candidate (validation) group and assumed to be untested.
      t2 <- setdiff(t2_hybrids, trn)

      # Remaining hybrids, which have either one or zero parents.
      t01 <- setdiff(x, union(t2, trn))

      # All hybrids, for which the Dent parent was an element of dt, but the
      # Flint parent was not an element of ft and vice versa, are assigned to
      # the T1 candidate (validation) group.
      t1_list <- hybrid_list[!t2_index]
      t1_index <- lapply(X = seq_len(length(t1_list)), FUN = function(j) {
        ind <- any(t1_list[[j]] %in% c(dent_trn, flint_trn))
        })
      t1_index <- unlist(t1_index)
      t1 <- t01[t1_index]

      # All hybrids in hybrid_vec, for which neither the Dent parent nor the
      # Flint parent was an element of dt or ft, respectively, were assigned to
      # the T0 candidate (validation) group.
      t0 <- setdiff(x, c(t2, t1, trn))
      n_t0 <- length(t0)
      n_t1 <- length(t1)
      n_t2 <- length(t2)
      n_trn <- length(trn)

      stopifnot(sum(n_t0, n_t1, n_t2, n_trn) == length(x))
    }
    cv <- utils::stack(list(TS = trn,
                            VS0 = t0,
                            VS1 = t1,
                            VS2 = t2))
    cv <- cv[order(cv$values), ]

    container[, i] <- as.character(cv$ind)
    if (i == 1) {
      rownames(container) <- as.character(cv$values)
      i <- i + 1
      # Only store a CV-run that is not a duplicate of previous runs.
    } else if (!any(duplicated(container[, 1:i], MARGIN = 2))) {
      i <- i + 1
    } else next
  }
 container
}

