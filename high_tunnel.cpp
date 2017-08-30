
// When fully debugged, this can be employed (it prevents bounds checks)
// #define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <sitmo.h>
#include <random>


using namespace Rcpp;
using namespace std;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo, sitmo)]]
//[[Rcpp::plugins(cpp11)]]

// 
// Expand from values per instar to per stage (i.e., day)
// 
//[[Rcpp::export]]
mat instar_to_stage(vec stage_values, uword n_stages, mat stage_days) {

    uword n_lines = stage_days.n_rows, n_instars = stage_days.n_cols;
    
    mat stage_days_cs = cumsum(stage_days, 1);
    
    mat out(n_lines, n_stages, fill::zeros);
    
    if (stage_values.n_elem != n_instars) {
        stop("ncol in stage_days should equal length of stage_values");
    }
    
    for (uword i = 0; i < n_lines; i++) {
        out(i, span(0, stage_days_cs(i,0)-1)).fill(stage_values[0]);
        for (uword j = 1; j < n_instars; j++) {
            out(i, span(stage_days_cs(i,j-1), stage_days_cs(i,j)-1)).fill(
                    stage_values[j]);
        }
    }
    
    return out;
}






// 
// Create Leslie matrix from aphid info
// 
//[[Rcpp::export]]
mat leslie_matrix(int n_stages, ivec stage_days, double surv_juv, 
                   vec surv_adult, vec repro) {
    vec tmp;
    mat LL;
    int juv_time = accu(stage_days(span(0, (stage_days.n_elem - 2))));
    // Age-specific survivals
    tmp = vec(n_stages - 1);
    tmp.head(juv_time).fill(surv_juv);
    tmp.tail(n_stages-juv_time-1) = surv_adult(span(0,(n_stages-juv_time-2)));
    LL = diagmat(tmp, -1);
    // Age-specific fecundities
    LL(0, span(juv_time, juv_time + stage_days(stage_days.n_elem - 1) - 1)) =
        repro(span(0, stage_days(stage_days.n_elem - 1) - 1)).t();
    
    return LL;
}



// This computes the "stable age distribution" from the Leslie matrix, which is 
// the proportion of different classes that is required for the population to grow 
// exponentially
//[[Rcpp::export]]
mat leslie_sad(mat L) {
    
    cx_vec r_cx;
    cx_mat SAD;
    
    eig_gen(r_cx, SAD, L);
    
    vec r = abs(r_cx);

    double rmax = max(r);
    
    cx_mat SADdist = SAD.cols(arma::find(r == rmax));
    cx_double all_SAD = accu(SADdist);
    SADdist /= all_SAD;
    SADdist.resize(SADdist.n_elem, 1);
    
    mat SADdist_Re = arma::real(SADdist);
    
    return SADdist_Re;
}



// Equation 6 from the paper
//[[Rcpp::export]]
mat attack_probs(double a, vec p_i, double Y_m, double x, double h, double k, 
                 vec resist_surv) {
    mat mm = (a * p_i * Y_m) / (h * x + 1);
    mat AA = (1 + mm / k);
    if (resist_surv.n_elem == 0) {
        AA = arma::pow(AA, -k);
    } else {
        if (resist_surv.n_elem != 2) stop("resist_surv should be of length 2");
        AA = arma::pow(AA, -k) + 
            resist_surv(0) * mm % arma::pow(AA, -k-1) + 
            resist_surv(1) * (1-(arma::pow(AA, -k) + mm % arma::pow(AA, -k-1)));
    }
    return AA;
}





// 
// Make new column of parasitoid abundances [i.e., Y(t+1), last 4 lines of Equation 2 
// in paper]; below, S_y_zt = S_y(z(t))
// In paper, sex_ratio = 1/2, pred_rate not present
// 
//[[Rcpp::export]]
vec parasitoid_abunds(double S_y_zt, vec A, mat L, vec X, vec Y_t, double s_i, 
                      double s_y, int m_1, double sex_ratio, double pred_rate) {

    mat tmp;
    // Y(t+1)
    vec Y_t1 = Y_t;
    // Y_1(t+1)
    tmp = S_y_zt * (1 - A).t() * (L * X);
    Y_t1(0) = tmp(0);  // <-- (0) is to let this know it's a double; it's already 1x1
    // Y_i(t+1) for (i = 1, ..., m_1)
    Y_t1(span(1, m_1)) = s_i * S_y_zt * Y_t.head(m_1);
    // Y_i(t+1) for (i = m_1+1, ..., m-1)
    Y_t1(span(m_1+1, Y_t1.n_elem-2)) = pred_rate * Y_t(span(m_1, Y_t.n_elem-3));
    // Y_m(t+1)
    Y_t1.tail(1) = s_y * Y_t.tail(1) + sex_ratio * Y_t(Y_t.n_elem-2);

    return Y_t1;
}



// 
// Dispersal (not in the paper)
// 
//[[Rcpp::export]]
mat dispersal(mat P_mat, double disp_rate, uvec disp_stages) {
    mat Pm = P_mat;
    double immigration;
    if (disp_stages.n_elem == 0) disp_stages = regspace<uvec>(0, Pm.n_rows - 1);
    for (auto i : disp_stages) {
        immigration = disp_rate * mean(Pm.row(i));
        Pm.row(i) *= (1 - disp_rate); // emigration
        Pm.row(i) += immigration;
    }
    return(Pm);
}