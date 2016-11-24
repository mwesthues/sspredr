#' Construction of ETA with complete predictor data
#'
#' \code{complete_eta} returns a list with BGLR model parameters.
#'
#' @param x A matrix with BLUES of the complete predictor.
#' @param geno NULL. Optional character vector with names corresponding to one
#'  of the two parental pools.
#' @param as_kernel logical. Should a feature matrix or a variance covariance
#'  matrix be returned?
#' @param is_pedigree logical. 'TRUE' if the input is a pedigree matrix.
#' @param bglr_model A character specifying the algorithm that shall be used \
#'  when calling \code{BGLR()}.
#'
#' @return A list with model parameters used as input in the \code{BGLR()} call.
#' @examples
#' data("hybrid_nms", "imp_snps")
#' geno <- vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 1,
#'                FUN.VALUE = character(1))
#' x <- imp_snps[rownames(imp_snps) %in% geno, ]
#' eta <- complete_eta(x = x, geno = geno, bglr_model = "BRR")
#' str(eta)
#' @importFrom magrittr %>%
#' @export
complete_eta <- function(x, geno, as_kernel = FALSE, is_pedigree = FALSE,
                         bglr_model) {
  # Input tests
  if (class(x) != "matrix") stop("'x' is not a matrix")
  if (class(geno) != "character") stop("'geno' is not a character vector")
  stopifnot(bglr_model %in% c("FIXED", "BRR", "BayesA", "BayesB", "BayesC",
                              "RKHS"))
  comgeno <- unique(geno)
  x <- x[rownames(x) %in% comgeno, ]
  # Sort
  x <- x[match(comgeno, rownames(x)), ]
  M <- x[, matrixStats::colVars(x) != 0]

  if (isTRUE(is_pedigree)) {
    M <- M[match(comgeno, rownames(M)), match(comgeno, colnames(M))]
    diag(M) <- diag(M) + 0.01
    L <- M %>% chol() %>% t()
  } else if (isTRUE(as_kernel)) {
    L <- build_kernel(M = M)
  } else {
    L <- M
  }

  # Design matrix to map GCA effects to corresponding parents in y.
  Z <- Matrix::sparse.model.matrix(~-1 + factor(geno),
                                   drop.unused.levels = FALSE)
  colnames(Z) <- gsub("factor\\(geno\\)", replacement = "", x = colnames(Z))
  Z <- Z[, match(comgeno, colnames(Z))]
  rownames(Z) <- geno

  MG <- as.matrix(Z %*% L)

  # Output tests
  stopifnot(identical(rownames(MG), geno))

  list(list(X = MG, model = bglr_model))
}
