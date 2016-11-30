#' Leave-one-out cross validation
#'
#' \code{run_loocv} returns a data frame with predicted values.
#'
#' @param Pheno A matrix with phenotypic BLUES.
#' @param ETA A list with ETA-objects for BGLR.
#' @param hybrid Boolean. Indicates whether prediction is done for hybrids.
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
#' eta <- impute_eta(x = x, y = y, as_kernel = TRUE, is_pedigree = FALSE,
#'                   geno = geno, bglr_model = "BRR")
#' eta[] <- lapply(seq_along(eta), FUN = function(i) {
#'   dat <- eta[[i]]
#'   x <- dat[["X"]]
#'   rownames(x) <- hybrid_nms
#'   dat[["X"]] <- x
#'   dat
#' })
#' dir.create("./bglr_out")
#' run_loocv(Pheno = Pheno, ETA = eta, hybrid = TRUE, father_idx = 2,
#'           mother_idx = 1, split_char = "_", trait = "gtm", iter = 2500L,
#'           speed_tst = FALSE, run = 1L, verbose = FALSE,
#'           out_loc = "./bglr_out/")
#' unlink("./bglr_out", recursive = TRUE)
#' @export
run_loocv <- function(Pheno, ETA, hybrid = TRUE,
                      father_idx = NULL, mother_idx = NULL,
                      split_char = NULL, trait, iter,
                      speed_tst = FALSE,
                      run = 1L, verbose = FALSE, out_loc) {

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
  stopifnot(typeof(iter) == "integer")
  stopifnot(typeof(speed_tst) == "logical")

  # Ensure that the elements of ETA and the phenotyipic values are in the same
  # order.
  eta_rownms <- lapply(ETA, FUN = function(x) {
    rownames(x[["X"]])
  })
  eta_rownms <- Reduce(intersect, eta_rownms)
  stopifnot(identical(rownames(Pheno), eta_rownms))

  # Hybrid progeny and parent genotypes names for LOOCV-matching.
  geno <- rownames(Pheno)
  if (isTRUE(hybrid)) {
    if (any(vapply(list(father_idx, mother_idx), FUN = is.null,
                   FUN.VALUE = logical(1)))) {
      stop("Specify position of paternal and maternal IDs for hybrids")
    }
    if (is.null(split_char)) {
      stop("Specify split character for hybrid names")
    }
    # Paternal IDs
    father <- vapply(strsplit(geno, split = split_char), FUN = "[[",
                     father_idx,
                     FUN.VALUE = character(1))
    # Maternal IDs
    mother <- vapply(strsplit(geno, split = split_char), FUN = "[[",
                     mother_idx,
                     FUN.VALUE = character(1))
    # Use only T0 hybrids for the training set.
    tst <- intersect(grep(mother[run], x = geno, invert = TRUE),
                     grep(father[run], x = geno, invert = TRUE))
  } else {
    tst <- grep(geno[run], x = geno)
  }
	y <- Pheno[, trait]
	y[-tst] <- NA_real_

  if (isTRUE(speed_tst)) {
    systime <- system.time(replicate(10, BGLR::BGLR(y = y,
                                                    ETA = ETA,
                                                    nIter = iter,
                                                    saveAt = out_loc,
                                                    burnIn = iter / 2,
                                                    verbose = verbose)))
    systime
  } else if (!isTRUE(speed_tst)) {

    # BGLR implementation
    res <- data.frame(Phenotype = NA_character_,
                      Geno = NA_character_,
                      Mother = NA_character_,
                      Father = NA_character_,
                      y = NA_complex_,
                      yhat = NA_complex_)

	  # run the model (GBLUP)
    mod_BGLR <- BGLR::BGLR(y = y,
                           ETA = ETA,
                           nIter = iter,
                           burnIn = iter / 2,
                           saveAt = out_loc,
                           verbose = verbose)

    # store results
    res$Phenotype <- trait
    res$Geno <- geno[run]
    if (isTRUE(hybrid)) {
      res$Mother <- mother[run]
      res$Father <- father[run]
    }
    res$y <- Pheno[run, trait]
    res$yhat <- mod_BGLR$yHat[run]
    res$var_yhat <- mod_BGLR$SD.yHat[run] ^ 2
    res
  }
}
