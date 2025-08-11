// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;

//' Two-sample Kolmogorov–Smirnov D statistic (fast)
 //' @noRd
 // [[Rcpp::export]]
 double ks2_fast(NumericVector x, NumericVector y) {
   int nx = x.size(), ny = y.size();
   std::sort(x.begin(), x.end());
   std::sort(y.begin(), y.end());
   int i = 0, j = 0;
   double cdf_x = 0.0, cdf_y = 0.0, dmax = 0.0;

   while (i < nx && j < ny) {
     if (x[i] <= y[j]) {
       ++i;   cdf_x = double(i) / nx;
     } else {
       ++j;   cdf_y = double(j) / ny;
     }
     double diff = std::fabs(cdf_x - cdf_y);
     if (diff > dmax) dmax = diff;
   }
   return dmax;
 }
