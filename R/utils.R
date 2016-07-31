#' Check TRN Assembly
#' @keywords internal
#' @param x A matrix containing cross-validation folds inside columns and
#'  hybrid names in rows.
#' @param hybrid_split Character that separates parents in hybrid names.
#' @return 'success', if all checks were successful.
check_cv <- function(x, hybrid_split = "_") {
  fl_nm <- names(x)
  trn_pos <- regexpr("(?<=_trn=)[0-9]+(?=_min_size=)", text = fl_nm, 
                     perl = TRUE)
  trn_nm <- substring(fl_nm, first = trn_pos, 
                      last = trn_pos + attr(trn_pos, "match.length") - 1)
  trn_nm <- as.numeric(trn_nm)
  min_pos <- regexpr("(?<=_min_size=)[0-9]+(?=_m=)", text = fl_nm, perl = TRUE)
  min_nm <- substring(fl_nm, first = min_pos,
                      last = min_pos + attr(min_pos, "match.length") - 1)
  min_nm <- as.numeric(min_nm)
  m_pos <- regexpr("(?<=_m=)[0-9]+(?=_f=)", text = fl_nm, perl = TRUE)
  m_nm <- substring(fl_nm, first = m_pos, 
                    last = m_pos + attr(m_pos, "match.length") - 1)
  m_nm <- as.numeric(m_nm)
  f_pos <- regexpr("(?<=_f=)[0-9]{2}", text = fl_nm, perl = TRUE)
  f_nm <- substring(fl_nm, first = f_pos, 
                    last = f_pos + attr(f_pos, "match.length") - 1)
  f_nm <- as.numeric(f_nm)
  dat <- x[[1]]

  unique_hybrids <- identical(length(unique(rownames(dat))),
                              nrow(dat))
  unique_runs <- !any(duplicated(dat, MARGIN = 2))
  if (!all(colSums(dat == "TRN") == trn_nm)) {
    stop("TRN size too small")
  }
  if (!all(colSums(dat == "T0") >= min_nm)) {
    stop("T0 size too small")
  }
  if (!all(colSums(dat == "T1") >= min_nm)) {
    stop("T1 size too small")
  }
  if (!all(colSums(dat == "T2") >= min_nm)) {
    stop("T2 size too small")
  }
  
  trn_parents <- vapply(seq_len(ncol(dat)), FUN = function(k) {
    trn_hyb <- rownames(dat)[dat[, k] == "TRN"]
    mother_par <- vapply(strsplit(trn_hyb, split = hybrid_split), 
                         FUN = "[[", 2, FUN.VALUE = character(1))
    father_par <- vapply(strsplit(trn_hyb, split = hybrid_split), 
                         FUN = "[[", 3, FUN.VALUE = character(1))
    all(c(all.equal(length(unique(mother_par)), m_nm),
          all.equal(length(unique(father_par)), f_nm)))
    }, FUN.VALUE = logical(1))
  if (!all(trn_parents)) {
    stop("Wrong parents")
  }

  trn_mat <- do.call(cbind, lapply(seq_len(ncol(dat)), FUN = function(j) {
    dat[, j] == "TRN"
  }))
  trn_duplicates <- sum(duplicated(trn_mat, MARGIN = 2))
  if (trn_duplicates != 0) stop("Duplicated TRN runs")

  t2_mat <- do.call(cbind, lapply(seq_len(ncol(dat)), FUN = function(j) {
    dat[, j] == "T2"
  }))
  t2_duplicates <- sum(duplicated(t2_mat, MARGIN = 2))
  if (t2_duplicates != 0) stop("Duplicated T2 runs")

  t1_mat <- do.call(cbind, lapply(seq_len(ncol(dat)), FUN = function(j) {
    dat[, j] == "T1"
  }))
  t1_duplicates <- sum(duplicated(t1_mat, MARGIN = 2))
  if (t1_duplicates != 0) stop("Duplicated T1 runs")

  t0_mat <- do.call(cbind, lapply(seq_len(ncol(dat)), FUN = function(j) {
    dat[, j] == "T0"
  }))
  t0_duplicates <- sum(duplicated(t0_mat, MARGIN = 2))
  if (t1_duplicates != 0) stop("Duplicated T0 runs")

  return("success")
}
