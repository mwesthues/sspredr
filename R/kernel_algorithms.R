#' Relationship Matrix
#'
#' \code{build_kernel} generates kinship matrices from raw feature matrices
#'
#' @param M A matrix. Genotype names are stored in rows whereas feature names
#'  are stored in columns.
#' @param lambda Numeric vector (scalar) to add to diagonal of the relationship
#'  matrix. This avoids numerical issues (singularities) when inverting the
#'  output.
#' @param algorithm Character vector (scalar) to generate the relationship
#'  matrix.
#' @param feat_weights Numeric vector specifying the call weights each feature
#'  has in the relationship matrix.
#' @return If \code{algorithm} is \code{RadenII} and \code{feat_weights} is set
#'  to \code{NULL}, features will be centered and scaled to unit variance;
#'  otherwise features will be centered and divided by
#'  \eqn{\sqrt{w_{i} * \sigma^{2}_{i}}}, where \eqn{w_{i}} is the weight of the
#'  \emph{i}th feature and \eqn{m_{i}} is the \emph{i}th feature. A typical
#'  use case for \emph{w} would be a vector of heritabilities for each feature.
#'
#'  If \code{algorithm} is \code{Zhang}, features will be centered and scaled
#'  by the sum of their variances.
#' @examples
#'  # Load a matrix with SNP genotypes encoded as numeric values
#'  data(mice)
#'  mice.snp <- mice.X[1:200, ]
#'
#'  # Generate a vector of feature weights.
#'  set.seed(90402)
#'  weights <- runif(n = ncol(mice.snp), min = 0, max = 1)
#'
#'  # Compute a SNP-based relationship matrix based on RadenII.
#'  radenII <- build_kernel(mice.snp, algorithm = "RadenII")
#'  dim(radenII)
#'  summary(c(radenII))
#'
#'  # Compute a SNP-based relationship matrix based on RadenII with weighted
#'  # features.
#'  weighted_radenII <- build_kernel(mice.snp, algorithm = "RadenII",
#'                                   feat_weights = weights)
#'
#'  #  Compute a SNP-based relationship matrix based on Zhang.
#'  zhang <- build_kernel(mice.snp, algorithm = "Zhang")
#'  all.equal(radenII, zhang)
#'  @export
build_kernel <- function(M, lambda = 0.01, algorithm = "RadenII",
                         feat_weights = NULL) {
 if (!typeof(M) %in% c("integer", "double")) {
    stop("The SNP data are neither of type 'integer' nor 'double'")
  }
  if (anyNA(M)) stop("NAs in M not allowed")
  if (algorithm == "RadenII") {
    if (!is.null(feat_weights)) {
      if (length(feat_weights) != ncol(M)) {
        stop(paste0("length of 'feat_weights' does not correspond to number of",
                    "features in 'M'"))
      }
      weighted <- sapply(seq_len(ncol(M)), FUN = function(i) {
        sqrt(feat_weights[i] * var(M[, i]))
      })
      W <- scale(M, center = TRUE, scale = weighted)
    } else {
      # scale with the standard deviation of original features
      W <- scale(M, center = TRUE, scale = TRUE)
    }
    # G = WW' / m, where m is the number of features
    G <- tcrossprod(W) / ncol(W)
  } else if (algorithm == "Zhang") {
    # center
    M_scaled <- scale(M, scale = FALSE)
    # get variances
    vars <- apply(M_scaled, MARGIN = 2, FUN = var)
    # compute crossproduct and add lambda
    G <- ((1 - lambda) * tcrossprod(M_scaled)) / sum(vars)
  }
  # Add a small number to the diagonal elements to avoid singularity in G.
  diag(G) <- diag(G) + lambda
  G
}
