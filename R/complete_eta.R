#' Construction of ETA with complete predictor data
#'
#' \code{complete_eta} returns a list with BGLR model parameters.
#'
#' @param x A matrix with BLUES of the complete predictor.
#' @param geno NULL. Optional character vector with names corresponding to one
#'  of the two parental pools.
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
#' @export
complete_eta <- function(x, geno, bglr_model) {
  # Input tests
  if (class(x) != "matrix") stop("'x' is not a matrix")
  if (class(geno) != "character") stop("'geno' is not a character vector")
  stopifnot(bglr_model %in% c("FIXED", "BRR", "BayesA", "BayesB", "BayesC",
                              "RKHS"))
  comgeno <- unique(geno)
  x <- x[rownames(x) %in% comgeno, ]
  # Sort
  x <- x[match(comgeno, rownames(x)), ]
  x <- x[, matrixStats::colVars(x) != 0]

  # Design matrix to map GCA effects to corresponding parents in y.
  Z <- Matrix::sparse.model.matrix(~-1 + factor(geno),
                                   drop.unused.levels = FALSE)
  colnames(Z) <- gsub("factor\\(geno\\)", replacement = "", x = colnames(Z))
  Z <- Z[, match(comgeno, colnames(Z))]
  rownames(Z) <- geno

  # x-kernel
  G <- build_kernel(x, lambda = 0.01, algorithm = "RadenII")
  MG <- as(Z %*% G, "matrix")

  # Output tests
  stopifnot(identical(colnames(MG), comgeno))
  stopifnot(identical(rownames(MG), geno))

  list(list(X = MG, model = bglr_model))
}
