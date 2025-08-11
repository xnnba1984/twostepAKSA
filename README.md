# twostepAKSA

This package implements the two-step test for predictive biomarkers with zero inflation introduced in the paper "A Two-Step Test for Zero-Inflated Biomarkers in Early-Phase Clinical Trials":
- **Part A (spike test)** compares treatment arms within the X=0 stratum.
- **Part B (tail test)** applies an averaged KS approach (AKSA) across biomarker prefixes among X>0 subjects.
- The two p-values are combined via **Fisher** or **Brown**; permutation calibration is used throughout.

## Installation

```r
install.packages("remotes")
remotes::install_github("xnnba1984/twostepAKSA")
```

## Minimal example

```r
library(twostepAKSA)
set.seed(1)

N  <- 90
pi0 <- 0.4
n0 <- round(N * pi0)
x  <- c(rep(0, n0), runif(N - n0))      # zero-inflated biomarker
x  <- sample(x)                          # shuffle
trt <- sample(rep(0:1, length.out=N))   # 1:1 randomization
y   <- rnorm(N)                          # standard normal outcome

# Spike-only effect for illustration
y[x == 0 & trt == 1] <- y[x == 0 & trt == 1] + 0.8

fit <- twostep_test(y, trt, x, nperm = 1000, combine = "brown", brown_nbrown = 200, seed = 123)
fit$pA; fit$pB; fit$p_combined
fit$method
```
