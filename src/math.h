#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace std;

arma::mat leslie_matrix(arma::uvec instar_days, double surv_juv,
                        arma::vec surv_adult, arma::vec repro);

arma::vec leslie_sad(arma::mat L);