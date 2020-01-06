#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace std;

// Logit and inverse logit
//[[Rcpp::export]]
arma::vec logit(arma::vec p) {
    arma::vec out = arma::log(p / (1-p));
    return out;
}
//[[Rcpp::export]]
arma::vec inv_logit(arma::vec a){
    arma::vec out = 1 / (1 + arma::exp(-a));
    return out;
}

//' Create Leslie matrix from aphid info
//' 
//' @param instar_days Integer vector of the number of stages (days) per aphid instar.
//' @param surv_juv Single numeric of daily juvenile survival.
//' @param surv_adult Numeric vector of aphid adult survivals by stage.
//' @param repro Numeric vector of aphid reproductive rates by stage.
//' 
//' 
//' @export
//[[Rcpp::export]]
arma::mat leslie_matrix(arma::uvec instar_days, double surv_juv,
                        arma::vec surv_adult, arma::vec repro) {
    uint n_stages = arma::sum(instar_days);
    arma::vec tmp;
    arma::mat LL;
    uint juv_time = arma::accu(instar_days(arma::span(0, (instar_days.n_elem - 2))));
    // Age-specific survivals
    tmp = arma::vec(n_stages - 1);
    tmp.head(juv_time).fill(surv_juv);
    tmp.tail(n_stages-juv_time-1) = surv_adult(arma::span(0,(n_stages-juv_time-2)));
    LL = arma::diagmat(tmp, -1);
    // Age-specific fecundities
    LL(0, arma::span(juv_time, juv_time + instar_days(instar_days.n_elem - 1) - 1)) =
        repro(arma::span(0, instar_days(instar_days.n_elem - 1) - 1)).t();
    
    return LL;
}


// This computes the "stable age distribution" from the Leslie matrix, which is
// the proportion of different classes that is required for the population to grow
// exponentially
// Used for constructing X_0 for a aphid_const class
arma::vec leslie_sad(arma::mat L) {
    
    arma::cx_vec r_cx;
    arma::cx_mat SAD;
    
    arma::eig_gen(r_cx, SAD, L);
    
    arma::vec r = arma::abs(r_cx);
    
    double rmax = arma::max(r);
    
    arma::cx_mat SADdist = SAD.cols(arma::find(r == rmax));
    arma::cx_double all_SAD = arma::accu(SADdist);
    SADdist /= all_SAD;
    SADdist.resize(SADdist.n_elem, 1);
    
    arma::vec SADdist_Re = arma::real(SADdist);
    
    return SADdist_Re;
}
