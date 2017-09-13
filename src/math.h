# ifndef _MATH_H
# define _MATH_H


#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace std;

arma::vec logit(arma::vec p);
arma::vec inv_logit(arma::vec a);

arma::mat leslie_matrix(arma::uvec instar_days, double surv_juv,
                        arma::vec surv_adult, arma::vec repro);

arma::vec leslie_sad(arma::mat L);

#endif