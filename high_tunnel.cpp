#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

//[[Rcpp::depends(RcppArmadillo)]]

// 
// Equivalent to diag(v, k) in MATLAB
// 
//[[Rcpp::export]]
arma::mat diag_mat(arma::vec v, int k = 0) {
    arma::mat dm = arma::diagmat(v, k);
    return dm;
}