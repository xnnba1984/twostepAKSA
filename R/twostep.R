#' @title Two-step test for predictive biomarkers with a zero spike
#' @description
#' Implements the two-step permutation test used in your simulations:
#' Part A targets the spike stratum (X = 0); Part B is the AKSA-style
#' tail scan over the positive biomarker values (X > 0). The two p-values
#' are combined by Fisher or Brown.
#'
#' @param y Numeric outcome vector (length N).
#' @param trt Integer or logical treatment indicator (0/1; length N).
#' @param x Numeric biomarker (length N), allowed to be zero-inflated.
#' @param nperm Integer. Number of random permutations for the reference
#'   distributions in Part A and Part B. Default 1000.
#' @param combine Character. "fisher" (default) or "brown".
#' @param brown_nbrown Integer. When `combine = "brown"`, the correlation
#'   adjustment is estimated by a light null-permutation scheme with this
#'   many replicates (default 200). Set to 0 to skip and default to Fisher.
#' @param min_pos Integer. Minimum number of positive-biomarker subjects
#'   required to run Part B (default 5). If fewer, Part B returns p=1.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list with elements \itemize{
#'   \item `pA` – Part A p-value (spike test).
#'   \item `pB` – Part B p-value (tail scan).
#'   \item `p_combined` – combined p-value by Fisher or Brown.
#'   \item `method` – combination method used.
#'   \item `details` – list with auxiliary information.
#' }
#' @export
twostep_test <- function(y, trt, x,
                         nperm = 1000,
                         combine = c("fisher","brown"),
                         brown_nbrown = 200,
                         min_pos = 5,
                         seed = NULL) {
  combine <- match.arg(combine)
  stopifnot(length(y) == length(trt), length(trt) == length(x))
  if (!is.null(seed)) set.seed(seed)

  y <- as.numeric(y)
  trt <- as.integer(trt)
  x <- as.numeric(x)

  # ----- Part A: spike test (X == 0), difference in means with permutation
  safe_diff <- function(a, b) {
    if (length(a) == 0 || length(b) == 0) return(0)
    abs(mean(a) - mean(b))
  }

  if (all(x > 0)) {
    pA <- 1
    statA <- 0
  } else {
    # Label G encodes (B==0,T) vs everything else; then compare groups 1 vs 0
    G <- trt + 2L * (x > 0)
    statA <- safe_diff(y[G == 1L], y[G == 0L])
    permA <- replicate(nperm, {
      gperm <- sample(G)
      safe_diff(y[gperm == 1L], y[gperm == 0L])
    })
    pA <- (1 + sum(permA >= statA)) / (nperm + 1)
  }

  # ----- Part B: tail scan (X > 0) averaged KS over sorted prefixes
  pos <- which(x > 0)
  if (length(pos) < min_pos) {
    pB <- 1
    statB <- 0
  } else {
    ord  <- order(x[pos])
    Ypos <- y[pos][ord]
    Tpos <- trt[pos][ord]
    npos <- length(Ypos)

    ks_vec <- vapply(1:(npos - 1), function(k) {
      yk <- Ypos[1:k]; tk <- Tpos[1:k]
      if (all(tk == 0) || all(tk == 1)) return(0)
      ks2_fast(yk[tk == 1], yk[tk == 0])
    }, numeric(1))
    statB <- mean(ks_vec)

    permB <- replicate(nperm, {
      tperm <- sample(Tpos)
      d <- vapply(1:(npos - 1), function(k) {
        yk <- Ypos[1:k]; tk <- tperm[1:k]
        if (all(tk == 0) || all(tk == 1)) return(0)
        ks2_fast(yk[tk == 1], yk[tk == 0])
      }, numeric(1))
      mean(d)
    })
    pB <- (1 + sum(permB >= statB)) / (nperm + 1)
  }

  # ----- Combine
  S_fisher <- -2 * (log(pA) + log(pB))
  pF <- stats::pchisq(S_fisher, df = 4, lower.tail = FALSE)

  p_combined <- pF
  method <- "Fisher"

  if (combine == "brown" && brown_nbrown > 0) {
    # Estimate correlation under null by light joint permutations:
    # permute treatment globally, recompute pA and pB for each permuted dataset
    # (each of those p-values is computed with *no extra* internal permutations;
    #  we approximate the null joint variability by reusing the observed stats).
    # For stability, we compute correlation on the scale -2*log(p).
    null_s <- replicate(brown_nbrown, {
      trtp <- sample(trt)
      # Recompute Part A stat under this permuted treatment, but *without*
      # inner permutations—use the same statistic as above, then approximate
      # the p via one-step plug-in using `permA`/`permB` empirical CDFs.
      # For simplicity and stability, we re-run the A/B p-values with fewer perms.
      pA_b <- .p_spike_once(y, trtp, x, nperm = ceiling(nperm/5))
      pB_b <- .p_tail_once (y, trtp, x, nperm = ceiling(nperm/5), min_pos = min_pos)
      c(-2 * log(pA_b), -2 * log(pB_b))
    })
    S1 <- null_s[1, ]; S2 <- null_s[2, ]
    rho <- suppressWarnings(stats::cor(S1, S2, method = "spearman", use = "complete.obs"))
    if (is.na(rho)) rho <- 0
    cfac <- 1 + rho
    df_b <- 4 / cfac
    pBROWN <- stats::pchisq(S_fisher / cfac, df = df_b, lower.tail = FALSE)
    p_combined <- pBROWN
    method <- sprintf("Brown (rho≈%.3f)", rho)
  }

  list(
    pA = pA,
    pB = pB,
    p_combined = p_combined,
    method = method,
    details = list(statA = statA, statB = statB, nperm = nperm,
                   combine = combine, brown_nbrown = brown_nbrown,
                   min_pos = min_pos)
  )
}

# internal helpers used by Brown's null correlation estimation
.p_spike_once <- function(y, trt, x, nperm) {
  safe_diff <- function(a, b) {
    if (length(a) == 0 || length(b) == 0) return(0)
    abs(mean(a) - mean(b))
  }
  if (all(x > 0)) return(1)
  G <- trt + 2L * (x > 0)
  obs <- safe_diff(y[G == 1L], y[G == 0L])
  perm <- replicate(nperm, {
    gperm <- sample(G)
    safe_diff(y[gperm == 1L], y[gperm == 0L])
  })
  (1 + sum(perm >= obs)) / (nperm + 1)
}

.p_tail_once <- function(y, trt, x, nperm, min_pos) {
  pos <- which(x > 0)
  if (length(pos) < min_pos) return(1)
  ord  <- order(x[pos])
  Ypos <- y[pos][ord];  Tpos <- trt[pos][ord]
  npos <- length(Ypos)

  ks_vec <- vapply(1:(npos - 1), function(k) {
    yk <- Ypos[1:k]; tk <- Tpos[1:k]
    if (all(tk == 0) || all(tk == 1)) return(0)
    ks2_fast(yk[tk == 1], yk[tk == 0])
  }, numeric(1))
  obs <- mean(ks_vec)

  permB <- replicate(nperm, {
    tperm <- sample(Tpos)
    d <- vapply(1:(npos - 1), function(k) {
      yk <- Ypos[1:k]; tk <- tperm[1:k]
      if (all(tk == 0) || all(tk == 1)) return(0)
      ks2_fast(yk[tk == 1], yk[tk == 0])
    }, numeric(1))
    mean(d)
  })
  (1 + sum(permB >= obs)) / (nperm + 1)
}
