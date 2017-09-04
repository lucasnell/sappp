
// When fully debugged, this can be employed (it prevents bounds checks)
// #define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <sitmo.h>
#include <cmath>
#include <random>


using namespace Rcpp;
using namespace std;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo, sitmo)]]
//[[Rcpp::plugins(cpp11)]]



// 
// Equation 6 from the paper
// Note: If attack_surv has length > 2, the 3rd+ items are ignored
// If attack_surv has length < 2, it's ignored entirely
// 
//[[Rcpp::export]]
vec attack_probs(double a, vec p_i, double Y_m, double x, double h, double k, 
                 vec attack_surv) {
    vec mm = (a * p_i * Y_m) / (h * x + 1);
    vec AA = (1 + mm / k);
    if (attack_surv.n_elem < 2 || sum(attack_surv) == 0) {
        AA = arma::pow(AA, -k);
    } else {
        AA = arma::pow(AA, -k) + 
            attack_surv(0) * mm % arma::pow(AA, -k-1) + 
            attack_surv(1) * (1-(arma::pow(AA, -k) + mm % arma::pow(AA, -k-1)));
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



// 
// Add process error to Y and X values
// Note: sigma_d_mult is a way to scale demographic stochasticity
// 
//[[Rcpp::export]]
List process_error(vec X_t1, vec Y_t1, vec X_t, vec Y_t, 
                   double sigma_x,  double sigma_y, double rho, 
                   double z, double Y_m, 
                   uword total_stages, uword living_aphids, 
                   double sigma_d_mult = 1) {
    
    mat Se(total_stages, total_stages, fill::zeros);
    
    if (sigma_d_mult == 0 && sigma_x == 0 && sigma_y == 0) {
        return List::create(Named("aphids") = X_t1, Named("wasps") = Y_t1);
    }
    
    // Aphid (both parasitized and not) process error
    Se(span(0, living_aphids-1), span(0, living_aphids-1)) =
        (sigma_x*sigma_x + sigma_d_mult * std::min(0.5, 1 / std::abs(1 + z))) *
        (rho * mat(living_aphids,living_aphids,fill::ones) + 
        (1-rho) * mat(living_aphids,living_aphids,fill::eye));
    // Version Tony sent (above is paper version):
    // sigma_x^2 * (rho * matrix(1,living_aphids,living_aphids) + (1-rho) * 
    // diag(1,living_aphids))
    
    // Mummy process error, turning back to zero
    Se(span(living_aphids,Se.n_rows-2), span(living_aphids,Se.n_cols-2)).fill(0);
    Se(span(0,living_aphids-1), span(living_aphids,Se.n_cols-2)).fill(0);
    Se(span(living_aphids,Se.n_rows-2), span(0,living_aphids-1)).fill(0);
    
    // Adult parasitoid process error
    Se(Se.n_rows-1, Se.n_cols-1) = sigma_y*sigma_y + 
        sigma_d_mult * std::min(0.5, 1 / std::abs(1 + Y_m));
    
    // chol doesn't work with zeros on diagonal
    uvec non_zero = find(Se.diag() > 0);
    
    // Cholesky decomposition of Se so output has correct variance-covariance matrix
    //   "a vector of independent normal random variables,
    //   when multiplied by the transpose of the Cholesky deposition of [Se] will
    //   have covariance matrix equal to [Se]."
    mat chol_decomp = chol(Se(non_zero,non_zero)).t();
    
    // Random numbers from distribution N(0,1)
    vec rnd = rnorm(non_zero.n_elem);
    
    // Making each element of rnd have correct variance-covariance matrix
    vec E = chol_decomp * rnd;
    
    uvec nz_aphid = non_zero(find(non_zero < X_t1.n_rows));
    uvec nz_wasp = non_zero(find(non_zero >= X_t1.n_rows)) - X_t1.n_rows;
    
    mat X_t1_e = X_t1;
    mat Y_t1_e = Y_t1;
    
    X_t1_e.rows(nz_aphid) = X_t1_e.rows(nz_aphid) % exp(E.head(nz_aphid.n_elem));
    Y_t1_e.rows(nz_wasp) = Y_t1_e.rows(nz_wasp) % exp(E.tail(nz_wasp.n_elem));
    
    // Because we used normal distributions to approximate demographic and environmental 
    // stochasticity, it is possible for aphids and parasitoids to 
    // "spontaneously appear" when the estimate of e(t) is large. To disallow this 
    // possibility, the number of aphids and parasitized aphids in a given age class 
    // on day t was not allowed to exceed the number in the preceding age class on 
    // day t – 1.
    
    for (uword i = 1; i < X_t1_e.n_elem; i++) {
        X_t1_e(i) = std::min(X_t1_e(i), X_t(i-1));
        if (X_t1_e(i) > X_t(i-1)) X_t1_e(i) = X_t(i-1);
    }
    // Not going to the end for parasitoids bc you can have more adults than mummies
    // bc adults stay in that stage for multiple days
    for (uword i = 1; i < (Y_t1_e.n_elem-1); i++) {
        Y_t1_e(i) = std::min(Y_t1_e(i), Y_t(i-1));
        if (Y_t1_e(i) > Y_t(i-1)) Y_t1_e(i) = Y_t(i-1);
    }
    
    return List::create(Named("aphids") = X_t1_e, Named("wasps") = Y_t1_e);
}

