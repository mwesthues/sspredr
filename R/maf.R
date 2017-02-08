#' Minor Allele Frequency
#'
#' \code{compute_maf} returns results from minor allele frequency checks
#'
#' @param x A matrix. Genotype names are stored in rows whereas marker names
#'  are stored in columns.
#' @param output Character vector with output options.
#' @param missing Character vector providing the encoding of missing elements.
#' @param maf_threshold Numeric or complex vector specifying the minor allele
#'  frequency threshold.
#' @return If \code{output} is "marker_names" a character vector with marker
#'  names that have passed the quality check will be returned. If \code{output}
#'  is "marker_maf" a numeric vector with the minor allele frequency at each
#'  marker locus will be returned. If \code{output} is "recoded" a matrix with
#'  the original dimensions will be returned where the most frequent genotype at
#'  each locus will be encoded as the highest value and the least frequent
#'  genotype will be encoded as the smalles value in \code{x}.
#' @examples
#'  # Load a matrix with SNP genotypes encoded as numeric values
#'  data(marker_numeric)
#'
#'  # Return the names of all polymorphic markers.
#'  compute_maf(marker_numeric, output = "marker_names", missing = NA_real_,
#'              mafThresh = 0)
#'
#'  # Set one locus to monomorphic
#'  marker_mono <- marker_numeric
#'  marker_mono[, 1] <- 1
#'  compute_maf(marker_numeric, output = "marker_names", missing = NA_real_,
#'              mafThresh = 0)
#'
#'  # Load a matrix with SNP genotypes encoded as character values
#'  data(marker_character)
#'
#'  # Return the minor allele frequency at each locus with MAF >= 0.05.
#'  compute_maf(marker_character, output = "marker_maf", missing = "??",
#'              mafThresh = 0.05)
#'
#'  # Return the genotype of the major as well as of the minor allele at each
#'  # locus.
#'  compute_maf(marker_character, output = "recoded", missing = "??",
#'              mafThresh = 0)
#' @export
compute_maf <- function(
  x,
  output = c("marker_names", "marker_maf", "recoded"),
  missing,
  maf_threshold) {


  # Store the original type of the input so that results for the
  # 'output'-option 'recoded' can be converted back to the original type.
  x_type <- typeof(x)

  # General input checks
  if (class(x) != "matrix") {
    stop("Input has to be of class 'matrix'")
  }

  if (typeof(missing) != typeof(x)) {
    stop("Types of 'missing' and 'x' must be the same.")
  }

  if (any(colMeans(is.na(x)) > 0.1)) {
  warning("Some markers have more than 10% missing values")
  }

  unique_elements <- sort(na.exclude(unique(c(x))))

  # Checks that are specific to numeric input matrices.
  if (typeof(x) %in% c("double", "integer")) {
    # Get the unique numeric elements in order to check whether the encoding of
    # genotypes with respect to homozygous (recessive and dominant) and
    # heterozygous genotypes makes sense.
    unique_elements <- sort(na.exclude(unique(c(x))))
    if (length(unique_elements) > 3) {
      stop("Only recessive, dominant and heterozygous genotypes are allowed.")
    }
    if (length(unique_elements) < 2) {
      stop("Only one unique genotype present")
    }
    if (length(unique_elements) == 3) {
      if (!unique_elements[2] == mean(c(unique_elements[1], unique_elements[3]))) {
        stop("Heterozygous value must be the mean of the homozygous values")
      }
    }
  }


  # It is assumed that the minor genotype in 'x' is encoded with a lower value
  # than the major genotype and that the heterozygous genotype is encoded as
  # the arithmetic mean of the two values.
  if (length(unique_elements) == 3) {
    min_value <- unique_elements[1]
    het_value <- unique_elements[2]
    max_value <- unique_elements[3]
    n_het <- colSums(x == het_value, na.rm = TRUE)
  } else {
    min_value <- unique_elements[1]
    max_value <- unique_elements[2]
  }

  # Count the occurrences of each genotype to determine the minor allele
  # frequency at each locus.
  n_min <- colSums(x == min_value, na.rm = TRUE)
  n_maj <- colSums(x == max_value, na.rm = TRUE)

  if(typeof(x) == "character") {
    # For each locus, check whether there are no more than two alleles.
    x[x == missing] <- NA_character_
    if (length(unique_elements) > 3) {
      stop("Only recessive, dominant and heterozygous genotypes are allowed.")
    }
    if (length(unique_elements) < 2) {
      stop("Only one unique genotype present")
    }
    if (length(unique_elements) == 3) {
      dom_allele <- purrr::flatten_chr(strsplit(max_value, split = ""))[1]
      rec_allele <- purrr::flatten_chr(strsplit(min_value, split = ""))[1]
      split_het_value <- purrr::flatten_chr(strsplit(het_value, split = ""))
      if (!dom_allele %in% split_het_value || !rec_allele %in% split_het_value) {
        stop("Heterozygotes must be encoded by the unique letters of homozygotes")
      }
    }
  }


  recoding_vector <- rep(min_value, times = ncol(x))
  recoding_vector[n_maj >= n_min] <- max_value
  storage.mode(recoding_vector) <- "character"
  storage.mode(x) <- "character"

  # For each locus, determine whether the assumed major genotype is actually
  # the most frequent one. If this is not the case (if) recode the genotypes
  # to match the expected values in 'unique_elements', otherwise (else) leave
  # them as before.
  recoded_x <- lapply(seq_len(ncol(x)), FUN = function(i) {
    x_col <- x[, i]
    actual_major <- recoding_vector[i]
    if (actual_major != as.character(max_value)) {
      x_col[x_col == actual_major] <- "major"
      x_col[x_col == max_value] <- "minor"
    } else {
      x_col[x_col == min_value] <- "minor"
      x_col[x_col == max_value] <- "major"
    }
    x_col
  })
  recoded_x <- matrix(unlist(recoded_x), ncol = ncol(x), byrow = FALSE)
  recoded_x[recoded_x == "minor"] <- min_value
  recoded_x[recoded_x == "major"] <- max_value
  if (x_type %in% c("double", "integer")) {
    storage.mode(recoded_x) <- "integer"
  }

  # Count the occurrences of each genotype to determine the minor allele
  # frequency at each locus.
  n_min <- colSums(recoded_x == min_value, na.rm = TRUE)
  n_maj <- colSums(recoded_x == max_value, na.rm = TRUE)
  if (length(unique_elements) == 3) {
    n_het <- colSums(recoded_x == het_value, na.rm = TRUE)
  } else {
    n_het <- 0
  }
  # Number of genotypes, which are not missing.
  n_geno <- n_maj + n_min + n_het
  # Frequency of the major allele at any locus.
  p_maj <- ((2 * n_maj) + n_het) / (2 * n_geno)
  # Frequency of the minor allele at any locus.
  p_min <- 1 - p_maj

  output <- match.arg(output)
  switch(output,
         marker_names = colnames(x[, p_min >= maf_threshold]),
         marker_maf = list(major = p_maj, minor = p_min),
         recoded = recoded_x)
}
