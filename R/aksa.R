#' @title Raw AKSA (full-biomarker scan) with permutation p-value
#' @description Sorts by the biomarker and averages two-sample KS distances
#' across prefixes (the original AKSA scan); permutes treatment to obtain
#' an exact reference under the null.
#' @inheritParams twostep_test
#' @return A list with fields `stat`, `pval`, and `details`.
#' @export
aksa_test <- function(y, trt, x, nperm = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ord  <- order(x)
  Y <- y[ord];  T <- trt[ord];  n <- length(Y)
  ks_vec <- vapply(1:(n - 1), function(k) {
    yk <- Y[1:k]; tk <- T[1:k]
    if (all(tk == 0) || all(tk == 1)) return(0)
    ks2_fast(yk[tk == 1], yk[tk == 0])
  }, numeric(1))
  stat <- mean(ks_vec)

  perm <- replicate(nperm, {
    tperm <- sample(T)
    d <- vapply(1:(n - 1), function(k) {
      yk <- Y[1:k]; tk <- tperm[1:k]
      if (all(tk == 0) || all(tk == 1)) return(0)
      ks2_fast(yk[tk == 1], yk[tk == 0])
    }, numeric(1))
    mean(d)
  })
  pval <- (1 + sum(perm >= stat)) / (nperm + 1)
  list(stat = stat, pval = pval, details = list(nperm = nperm))
}
