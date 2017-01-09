utils::globalVariables(c("."))
#' Imputation of mRNAs via pedigree and genomic information.
#'
#' \code{impute2} returns a list with BGLR model parameters.
#'
#' @param ped A matrix with pedigree records.
#' @param snp A matrix with SNP records.
#' @param mrna A matrix with transcriptomic BLUEs.
#' @param as_kernel logical. Should a feature matrix or a variance covariance
#'  matrix be returned? Should be 'FALSE' if 'x' denotes pedigree information.
#' @param geno NULL. Optional character vector with names corresponding to one
#'  of the two parental pools.
#' @param bglr_model A character specifying the algorithm that shall be used \
#'  when calling \code{BGLR()}.
#'
#' @return A list with model parameters used as input in the \code{BGLR()} call.
#' @examples
#' data("mice_ped", "mice_snp", "mice_mrna")
#' geno <- rownames(mice_ped)
#' eta <- impute2(ped = mice_ped, snp = mice_snp, mrna = mice_mrna, geno = geno,
#'                bglr_model = "BRR")
#' str(eta)
#' @importFrom magrittr %>%
#' @export
impute2 <- function(ped, snp, mrna, as_kernel = TRUE, geno, bglr_model) {
  # Input tests
  stopifnot(isTRUE(is.matrix(ped)))
  stopifnot(isTRUE(is.matrix(snp)))
  stopifnot(isTRUE(is.matrix(mrna)))

  # Specify which genotypes belong to pedigree (0), SNP (1) and mRNA (2)
  # records.
  category_tb <- tibble::tibble(G = unique(geno),
                                Category = 0)
  category_tb[category_tb$G %in% rownames(snp), "Category"] <- 1
  category_tb[category_tb$G %in% rownames(mrna), "Category"] <- 2
  category_df <- as.data.frame(category_tb)
  nm0 <- category_df[category_df$Category == 0, "G"]
  nm1 <- category_df[category_df$Category == 1, "G"]
  nm2 <- category_df[category_df$Category == 2, "G"]

  pred_lst <- list(ped = ped, snp = snp, mrna = mrna) %>%
    purrr::map(function(x) x[rownames(x) %in% geno, ]) %>%
    purrr::map_at(.at = "ped",
                  .f = function(x) x[, colnames(x) %in% geno])


  # Z ---------------------------------------------------------------------
  # Create design matrices associating the random effects with the response
  # variable.
  geno_lst <- tibble::tibble(G = geno) %>%
    dplyr::left_join(category_tb, by = "G") %>%
    split(.$Category) %>%
    purrr::map("G")
  names(geno_lst) <- c("ped", "snp", "mrna")
  geno_lvls <- geno_lst %>%
    purrr::map(factor) %>%
    purrr::map(levels)
  Z_lst <- geno_lst %>%
    purrr::map(~Matrix::sparse.model.matrix(~-1 + factor(.),
                                            drop.unused.levels = FALSE)) %>%
    purrr::map2(.y = geno_lst,
                .f = function(x, y) {
      rownames(x) <- y
      x
    }) %>%
    purrr::map2(.y = geno_lvls,
                .f = function(x, y) {
      colnames(x) <- y
      x
    })


  # J & M -----------------------------------------------------------------
  # For numerical reasons, add a small number to the diagonal of the pedigree
  # matrix. This allows the Cholesky decomposition to work.
  A <- pred_lst$ped %>%
    .[match(c(nm0, nm1, nm2), rownames(.)),
      match(c(nm0, nm1, nm2), colnames(.))]
  # Numerator relationship matrix.
  diag(A) <- diag(A) + 0.01
  Ainv <- solve(A)

  # Genomic relationship matrix.
  G <- pred_lst$snp %>%
    .[match(c(nm1, nm2), rownames(.)), ] %>%
    .[, matrixStats::colVars(.) != 0] %>%
    scale() %>%
    (function(x) tcrossprod(x) / ncol(x))
  diag(G) <- diag(G) + 0.01
  Ginv <- solve(G)

  # Quality check for transcriptomic data.
  M2 <- pred_lst$mrna %>%
    .[match(nm2, rownames(.)), ] %>%
    .[, matrixStats::colVars(.) != 0]

  A02 <- A[match(nm0, rownames(A)), match(nm2, colnames(A))]
  A22 <- A[match(nm2, rownames(A)), match(nm2, colnames(A))]
  M0 <- A02 %*% solve(A22) %*% M2

  G12 <- G[match(nm1, rownames(G)), match(nm2, colnames(G))]
  G22 <- G[match(nm2, rownames(G)), match(nm2, colnames(G))]
  M1 <- G12 %*% solve(G22) %*% M2

  J2 <- matrix(-1, nrow = length(nm2), ncol = 1)
  J1 <- G12 %*% solve(G22) %*% J2
  J0 <- A02 %*% solve(A22) %*% J2
  J_lst <- list(ped = J0, snp = J1, mrna = J2)


  # W ----------------------------------------------------------------------
  W <- Z_lst %>%
    purrr::map2(.y = list(ped = M0, snp = M1, mrna = M2),
         .f = function(x, y) {x %*% y}) %>%
    purrr::map(as.matrix) %>%
    do.call(rbind, .)
  if (isTRUE(as_kernel)) {
    W <- W %>% .[match(c(nm0, nm1, nm2), rownames(.)), ] %>%
      build_kernel(lambda = 0.01, algorithm = "RadenII")
  }
  W <- W[match(geno, rownames(W)), ]

  # X* --------------------------------------------------------------------
  X_prime <- purrr::map2(.x = Z_lst, .y = J_lst, .f = function(x, y) x %*% y) %>%
    do.call(rbind, .)


  # U0 --------------------------------------------------------------------
  epsilon0 <- A %>%
    solve() %>%
    .[match(nm0, rownames(.)), match(nm0, colnames(.))] %>%
    solve() %>%
    chol() %>%
    t()
  U0_top <- Z_lst$ped %*% epsilon0
  U0_bottom <- matrix(0,
                      nrow = geno_lst %>%
                        .[c("snp", "mrna")] %>%
                        purrr::flatten() %>%
                        purrr::as_vector() %>%
                        length(),
                      ncol = ncol(U0_top))
  U0 <- rbind(U0_top,
              U0_bottom)
  rownames(U0) <- geno_lst %>%
    purrr::flatten() %>%
    purrr::as_vector()


  # U1 --------------------------------------------------------------------
  epsilon1 <- G %>%
    solve() %>%
    .[match(nm1, rownames(.)), match(nm1, colnames(.))] %>%
    solve() %>%
    chol() %>%
    t()
  U1 <- Z_lst$snp %*% epsilon1
  U1_top <- matrix(0,
                   nrow = geno_lst %>% .[["ped"]] %>% length(),
                   ncol = ncol(U1))
  U1_bottom <- matrix(0,
                      nrow = geno_lst %>% .[["mrna"]] %>% length(),
                      ncol = ncol(U1))
  U1 <- do.call(rbind, list(U1_top, U1, U1_bottom))
  rownames(U1) <- geno_lst %>%
    purrr::flatten() %>%
    purrr::as_vector()


  # Output tests ----------------------------------------------------------
  # Make sure that all ETA-entries are sorted correctly.
  param_lst <- list(X = X_prime, W = W, U0 = U0, U1 = U1)
  param_lst <- param_lst %>%
    purrr::map(function(x) {
      x[match(geno, rownames(x)), ]
    })
  stopifnot(1 ==
    param_lst %>%
      purrr::map(rownames) %>%
      unique() %>%
      length()
  )


  # Return the ETA objects.
  list(list(X = param_lst$X, model = "FIXED"),
       list(X = param_lst$W, model = bglr_model),
       list(X = param_lst$U0, model = bglr_model),
       list(X = param_lst$U1, model = bglr_model))
}
