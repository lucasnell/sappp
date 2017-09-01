#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace std;

//
// Create Leslie matrix from aphid info
//
arma::mat leslie_matrix(arma::uvec stage_days, double surv_juv,
                        arma::vec surv_adult, arma::vec repro) {
    uint n_stages = arma::sum(stage_days);
    arma::vec tmp;
    arma::mat LL;
    uint juv_time = arma::accu(stage_days(arma::span(0, (stage_days.n_elem - 2))));
    // Age-specific survivals
    tmp = arma::vec(n_stages - 1);
    tmp.head(juv_time).fill(surv_juv);
    tmp.tail(n_stages-juv_time-1) = surv_adult(arma::span(0,(n_stages-juv_time-2)));
    LL = arma::diagmat(tmp, -1);
    // Age-specific fecundities
    LL(0, arma::span(juv_time, juv_time + stage_days(stage_days.n_elem - 1) - 1)) =
        repro(arma::span(0, stage_days(stage_days.n_elem - 1) - 1)).t();
    
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
