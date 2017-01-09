utils::globalVariables(c("."))
#' Imputation of ETA
#'
#' \code{impute_eta} returns a list with BGLR model parameters.
#'
#' @param x A matrix with BLUES of the complete predictor.
#' @param y A matrix with BLUEs of the incomplete predictor.
#' @param geno NULL. Optional character vector with names corresponding to one
#'  of the two parental pools.
#' @param as_kernel logical. Should a feature matrix or a variance covariance
#'  matrix be returned?
#' @param is_pedigree logical. Should be set to 'TRUE' if 'x' is a pedigree
#'  matrix.
#' @param bglr_model A character specifying the algorithm that shall be used \
#'  when calling \code{BGLR()}.
#'
#' @return A list with model parameters used as input in the \code{BGLR()} call.
#' @examples
#' data("hybrid_nms", "mrna", "imp_snps")
#' geno <- vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 1,
#'                FUN.VALUE = character(1))
#' x <- imp_snps[rownames(imp_snps) %in% geno, ]
#' y <- mrna[rownames(mrna) %in% geno, ]
#' eta <- impute_eta(x = x, y = y, as_kernel = TRUE, is_pedigree = FALSE, 
#'                   geno = geno, bglr_model = "BRR")
#' str(eta)
#' @importFrom magrittr %>%
#' @export
impute_eta <- function(x, y, as_kernel = TRUE, is_pedigree = TRUE, geno = NULL,
                       bglr_model) {
  # Input tests
  stopifnot(class(x) == "matrix")
  stopifnot(class(y) == "matrix")
  stopifnot(all(rownames(y) %in% rownames(x)))
  stopifnot(bglr_model %in% c("FIXED", "BRR", "BayesA", "BayesB", "BayesC",
                              "RKHS"))
  if (is.null(geno)) {
    geno <- rownames(x)
  }

  # Names of genotypes for which transcriptomic records exist
  nm2 <- rownames(y)
  # Names of genotypes for which transcriptomic records are missing.
  nm1 <- setdiff(rownames(x), nm2)

  # Hybrid parents not in y.
  geno1 <- geno[geno %in% nm1]
  # Hybrid parents in y.
  geno2 <- geno[geno %in% nm2]
  x <- x[match(unique(geno), rownames(x)), ]

  # Design matrix mapping GCA effects from genotypes without transriptomic
  # records to y1.
  Z1 <- Matrix::sparse.model.matrix(~-1 + factor(geno1),
                                    drop.unused.levels = FALSE)
  colnames(Z1) <- gsub("factor\\(geno1\\)", replacement = "", x = colnames(Z1))
  Z1 <- Z1[, match(nm1, colnames(Z1))]
  rownames(Z1) <- geno1
  # Design matrix mapping GCA effects from genotypes with transcriptomic
  # records to y2.
  Z2 <- Matrix::sparse.model.matrix(~-1 + factor(geno2),
                                    drop.unused.levels = FALSE)
  colnames(Z2) <- gsub("factor\\(geno2\\)", replacement = "", x = colnames(Z2))
  Z2 <- Z2[, match(nm2, colnames(Z2))]
  rownames(Z2) <- geno2

  y <- y[nm2, ]
  M2 <- y[, matrixStats::colVars(y) != 0]
  x <- x[match(c(nm1, nm2), rownames(x)),
         matrixStats::colVars(x) != 0]

  if (isTRUE(is_pedigree)) {
    A <- x[match(unique(geno), rownames(x)),
           match(unique(geno), colnames(x))]
    diag(A) <- diag(A) + 0.01
  } else {
    A <- tcrossprod(scale(x)) / ncol(x)
    diag(A) <- diag(A) + 0.01
  }
  A11 <- A[match(nm1, rownames(A)), match(nm1, colnames(A))]
  A12 <- A[match(nm1, rownames(A)), match(nm2, colnames(A))]
  A21 <- A[match(nm2, rownames(A)), match(nm1, colnames(A))]
  A22 <- A[match(nm2, rownames(A)), match(nm2, colnames(A))]
  Ainv <- solve(A)
  dimnames(Ainv) <- dimnames(A)
  A_up11 <- Ainv[match(nm1, rownames(Ainv)), match(nm1, colnames(Ainv))]
  A_up12 <- Ainv[match(nm1, rownames(Ainv)), match(nm2, colnames(Ainv))]
  # Eq.21
  M1 <- A12 %*% solve(A22) %*% M2
  stopifnot(identical(rownames(M1), nm1))
  J2 <- matrix(-1, nrow = ncol(A12), ncol = 1)
  # Eq.22
  J1 <- A12 %*% solve(A22) %*% J2
  # Eq.10
  epsilon <- t(chol(solve(A_up11)))

  # Eq.20
  W1 <- Z1 %*% M1
  W2 <- Z2 %*% M2
  W <- as.matrix(rbind(W1, W2))
  if (isTRUE(as_kernel)) {
    W <- W %>%
      .[match(c(nm1, nm2), rownames(.)), ] %>%
      build_kernel(lambda = 0.01, algorithm = "RadenII")
  }
  W <- W[match(geno, rownames(W)), ]
  U1 <- Z1 %*% epsilon
  U2 <- Z2 %*% matrix(0, nrow = ncol(Z2), ncol = ncol(U1))
  U <- as.matrix(rbind(U1, U2))
  U <- U[match(geno, rownames(U)), ]
  X_prime <- as.matrix(rbind(Z1 %*% J1, Z2 %*% J2))
  param_lst <- list(X = X_prime, W = W, U = U)

  stopifnot(all(unlist(lapply(param_lst, FUN = function(XWU) {
    identical(rownames(XWU), geno)
  }))))
  stopifnot(isTRUE(all(U2 == 0)))

  list(list(X = X_prime, model = "FIXED"),
       list(X = W, model = bglr_model),
       list(X = U, model = bglr_model))
}
