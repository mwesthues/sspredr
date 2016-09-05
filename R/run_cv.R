#' Cross validated prediction
#'
#' \code{run_cv} returns a data frame with predicted values.
#'
#' @param Pheno A matrix with phenotypic BLUES.
#' @param ETA A list with ETA-objects for BGLR.
#' @param cv A data.frame with the cross validation scheme.
#' @param father_idx NULL. A numeric or integer indicating the position of the
#'  parental ID in \code{hybrid}.
#' @param mother_idx NULL. A numeric or integer indicating the position of the
#'  maternal ID in \code{hybrid}.
#' @param split_char Character, which splits paternal and maternals IDs in
#'  \code{hybrid}.
#' @param trait Character vector with the name of the trait in \code{Pheno}.
#' @param iter Integer with the number of BGLR iterations.
#' @param speed_tst Boolean. Indicates whether a speed test (TRUE) or actual
#'  predictions (FALSE) shall be run.
#' @param run Integer specifying the current run corresponding to the predicted
#'  genotype in \code{Pheno}.
#' @param verbose Boolean. Shall BGLR output be printed to the console?
#' @param out_loc Path to the directory where temporary BGLR-output will be
#'  stored.
#'
#' @return Predicted values based on leave-one-out cross validation (LOOCV).
#' @examples
#' data("hybrid_nms", "mrna", "imp_snps", "Pheno")
#' geno <- vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 1,
#'                FUN.VALUE = character(1))
#' x <- imp_snps[rownames(imp_snps) %in% geno, ]
#' y <- mrna[rownames(mrna) %in% geno, ]
#' # BGLR ETA objects
#' eta <- impute_eta(x = x, y = y, geno = geno, bglr_model = "BRR")
#' eta[] <- lapply(seq_along(eta), FUN = function(i) {
#'   dat <- eta[[i]]
#'   x <- dat[["X"]]
#'   rownames(x) <- hybrid_nms
#'   dat[["X"]] <- x
#'   dat
#' })
#' # CV scheme generation
#' cv_mat <- sample_cv(hybrid_nms, n_mother = 39, n_father = 33,
#'                     n_hyb_trn = 200, min_size = 20, rounds = 2L,
#'                     hybrid_split = "_")
#' cv_mat <- cv_mat[[1]]
#' colnames(cv_mat) <- seq_len(ncol(cv_mat))
#' rownames(cv_mat) <- gsub("DF_", replacement = "", x = rownames(cv_mat))
#' cv <- as.data.frame(cv_mat)
#' cv$Sample_ID <- rownames(cv)
#' cv <- tidyr::gather(cv, Run, Set, -Sample_ID)
#' dir.create("./bglr_out")
#' # Run the prediction
#' run_cv(Pheno = Pheno, ETA = eta, cv = cv, father_idx = 2, mother_idx = 1,
#'        split_char = "_", trait = "gtm", iter = 2500L,
#'        speed_tst = FALSE, run = 1L, verbose = FALSE,
#'        out_loc = "./bglr_out/")
#' unlink("./bglr_out", recursive = TRUE)
#' @export
run_cv <- function(Pheno, ETA, cv, father_idx, mother_idx, split_char,
                   trait, iter, speed_tst, run, verbose = FALSE,
                   out_loc) {

  # Operations required for input checks.
  nrow_eta <- unique(vapply(ETA, FUN = function(x) {
    nrow(x[["X"]])
  }, FUN.VALUE = integer(1)))

  # Input testing.
  stopifnot(class(Pheno) == "matrix")
  stopifnot(class(ETA) == "list")
  stopifnot(length(nrow_eta) == 1)
  stopifnot(identical(nrow_eta, nrow(Pheno)))
  stopifnot(class(trait) == "character")
  stopifnot(length(trait) == 1)
  stopifnot(typeof(iter) %in% c("numeric", "integer"))
  stopifnot(typeof(speed_tst) == "logical")

  # Hybrid progeny and parent genotypes names for LOOCV-matching.
  hybrid <- rownames(Pheno)
  father <- vapply(strsplit(hybrid, split = split_char), FUN = "[[",
                   father_idx, FUN.VALUE = character(1))
  mother <- vapply(strsplit(hybrid, split = split_char), FUN = "[[",
                   mother_idx, FUN.VALUE = character(1))

  # Ensure that both, the phenotypic data and the CV sampling scheme are
  # complete.
  stopifnot(length(c(setdiff(cv$Sample_ID, hybrid),
                     setdiff(hybrid, cv$Sample_ID))) == 0)

  # Ensure that the elements of ETA and the phenotyipic values are in the same
  # order.
  eta_rownms <- lapply(ETA, FUN = function(x) {
    rownames(x[["X"]])
  })
  eta_rownms <- Reduce(intersect, eta_rownms)
  stopifnot(identical(rownames(Pheno), eta_rownms))
  cv_cur <- cv[cv$Run == run, ]
  cv_cur <- cv_cur[match(hybrid, cv_cur$Sample_ID), ]
  stopifnot(identical(rownames(Pheno), cv_cur$Sample_ID))

  # Get the indices of all test-set (T0, T1, T2) hybrids, define another vector
  # of thenotypic records and set the values of the latter to 'NA' if they belong
  # to a genotype that is part of the test set 'tst'.
  stopifnot(all(cv_cur$Set %in% c("T0", "T1", "T2", "TRN")))
  tst0 <- cv_cur$Set == "T0"
  tst1 <- cv_cur$Set == "T1"
  tst2 <- cv_cur$Set == "T2"
  tst <- as.logical(tst0 + tst1 + tst2)
  y <- yNA <- Pheno[, trait]
  yNA[tst] <- NA_real_

  if (isTRUE(speed_tst)) {
    systime <- system.time(replicate(10, BGLR::BGLR(y = yNA,
                                                    ETA = ETA,
                                                    nIter = iter,
                                                    burnIn = iter / 2,
                                                    saveAt = out_loc,
                                                    verbose = verbose)))
    stopifnot(class(systime) == "proc_time")
    systime
  } else if (!isTRUE(speed_tst)) {
    # RUN BGLR
    mod_BGLR <- BGLR::BGLR(y = yNA,
                           ETA = ETA,
                           nIter = iter,
                           burnIn = iter / 2,
                           saveAt = out_loc,
                           verbose = verbose)
    data.frame(Run = run,
               y = y,
               yHat = mod_BGLR$yHat,
               Set = cv_cur$Set)
  }
}
