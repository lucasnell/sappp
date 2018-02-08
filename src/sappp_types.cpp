
#include <RcppArmadillo.h> // arma namespace
#include <sitmo.h>         // parallel rng
#include <vector>          // vector class
#include <cmath>           // log, exp
#include <random>          // normal distribution
#include <cstdint>         // integer types
#include <algorithm>       // find
#include "math.h"          // leslie_matrix and leslie_sad
#include "sappp_types.h"     // all these classes


using namespace Rcpp;
using namespace std;

#define MAX_INT 4294967295



// Summarize a simulation's populations for all lines across all times and patches
// Rows are time, columns are lines, slices (3rd dimension) are patches

/* Overloaded fill function: */
// For one time, one line, one patch
void SimSummary::fill(uint t_, uint l_, uint p_, 
          double aphids_, double parasit_, double mummies_, double wasps_) {
    aphids(t_, l_, p_) = aphids_;
    parasit(t_, l_, p_) = parasit_;
    mummies(t_, l_, p_) = mummies_;
    wasps(t_, l_, p_) = wasps_;
}

// one time, one patch, ALL lines
void SimSummary::fill(uint t_, uint p_, 
          arma::rowvec aphids_, arma::rowvec parasit_, 
          arma::rowvec mummies_, arma::rowvec wasps_) {
    aphids.slice(p_).row(t_) = aphids_;
    parasit.slice(p_).row(t_) = parasit_;
    mummies.slice(p_).row(t_) = mummies_;
    wasps.slice(p_).row(t_) = wasps_;
}

// one time, ALL patches, ALL lines
// For the input matrices, rows are patches, columns are lines
void SimSummary::fill(uint t_,
          arma::mat aphids_, arma::mat parasit_,
          arma::mat mummies_, arma::mat wasps_) {
    for (uint i = 0; i < aphids.n_slices; i++) {
        aphids.slice(i).row(t_) = aphids_.row(i);
        parasit.slice(i).row(t_) = parasit_.row(i);
        mummies.slice(i).row(t_) = mummies_.row(i);
        wasps.slice(i).row(t_) = wasps_.row(i);
    }
}

arma::mat SimSummary::flatten() {
    
    uint n_patches = aphids.n_slices;
    uint n_lines = aphids.n_cols;
    uint n_t = aphids.n_rows;
    
    arma::mat output(n_lines * n_patches * n_t, 7);
    uint i = 0;
    
    for (uint p = 0; p < n_patches; p++) {
        for (uint l = 0; l < n_lines; l++) {
            for (uint t = 0; t < n_t; t++) {
                output(i,0) = static_cast<double>(p);
                output(i,1) = static_cast<double>(l);
                output(i,2) = static_cast<double>(t);
                output(i,3) = aphids(t,l,p);
                output(i,4) = parasit(t,l,p);
                output(i,5) = mummies(t,l,p);
                output(i,6) = wasps(t,l,p);
                i++;
            }
        }
    }
    
    return output;
}


void SimSummary::show() {
    
    uint n_t = aphids.n_rows;
    uint n_pops = aphids.n_cols;
    uint n_patches = aphids.n_slices;
    
    Rcout.precision(4);
    Rcout << std::fixed;
    
    Rcout << "< Info for SimSummary >" << endl;
    Rcout << "  time:    " << n_t       << endl;
    Rcout << "  pops:    " << n_pops    << endl;
    Rcout << "  patches: " << n_patches << endl;
    Rcout << endl;
    
    Rcout << "Fields:" << endl;
    Rcout << "  aphids:  " <<  "density of unparasitized aphids" << endl;
    Rcout << "  parasit: " <<  "density of parasitized, but alive, aphids" << endl;
    Rcout << "  mummies: " <<  "density of mummies" << endl;
    Rcout << "  wasps:   " <<  "density of adult wasps" << endl;
    Rcout << endl;
    
    Rcout << "Each field is a 3D array." << endl;
    Rcout << "Access elements like... " << endl;
    Rcout << "<SimSummary obj>$<field name>[<time>,<pop>,<patch>]" << endl;
}




// ======================================================================================
// ======================================================================================
// ======================================================================================

//          One aphid line and associated wasps

// ======================================================================================
// ======================================================================================
// ======================================================================================

// Using X and Y at time t+1 bc these get calculated before X and Y get iterated
// Aphid population
// Total (non-parasitized) aphids
double AphidPop::total_aphids() {
    return arma::sum(X_t1);
}

// Wasp population
// Total living, but parasitized aphids
double WaspPop::total_living_paras() {
    return arma::sum(Y_t1(arma::span(0, mum_days(0)-1)));
}
// Total adult wasps
double WaspPop::total_adult_wasps() {
    return Y_t1(Y_t1.n_elem-1);
}




// Wasp attack
// 
// Update attack probabilities
// Equation 6 from Meisner et al. (2014)
// Note: rel_attack is equivalent to p_i
// 
void WaspAttack::iterate_A(double Y_m, double x) {
    
    arma::vec mm = (a * rel_attack * Y_m) / (h * x + 1);
    arma::vec AA = (1 + mm / k);
    if (attack_surv.n_elem < 2 || sum(attack_surv) == 0) {
        A = arma::pow(AA, -k);
    } else {
        A = arma::pow(AA, -k) + 
            attack_surv(0) * mm % arma::pow(AA, -k-1) + 
            attack_surv(1) * (1-(arma::pow(AA, -k) + mm % arma::pow(AA, -k-1)));
    }
    return;
}






// Info about one aphid line and wasps that parasitize them
void AphidWasp::set_seed(uint seed) {
    rnd_engine.seed(seed);
}
// Harvest a patch (emulate harvesting alfalfa every month)
void AphidWasp::harvest() {
    // Kill non-parasitized aphids
    X_t1 *= harvest_surv;
    // Kill parasitized (but still living) aphids at the same rate
    Y_t1.head(mum_days(0)) *= harvest_surv;
    // Kill all mummies
    Y_t1(arma::span(mum_days(0), Y_t1.n_elem-2)).fill(0);
}

// Show an AphidWasp object info
void AphidWasp::show() const {
    arma::rowvec survs = arma::diagvec(leslie, -1).t();
    arma::rowvec fecs = leslie(0, arma::span(1, leslie.n_cols - 1));
    
    arma::uword N = std::min(static_cast<arma::uword>(6), survs.n_elem);
    
    Rcout.precision(4);
    Rcout << std::fixed;
    
    Rcout << "< Info for '" << aphid_name << "' aphid-wasp populations>" << endl;
    
    Rcout << "Aphid-line info (static):" << endl;
    Rcout << "  - resistances: (";
    Rcout << attack_surv(0) << ' ' << attack_surv(1) << ')' << endl;
    
    Rcout << "  - survivals:   (";
    for (unsigned i = 0; i < (N-1); i++) Rcout << survs(i) << ' ';
    if (survs.n_elem > N) Rcout << "... ";
    Rcout << survs(survs.n_elem-1) << ')' << endl;
    
    Rcout << "  - fecundities: (";
    for (unsigned i = 0; i < (N-1); i++) Rcout << fecs(i) << ' ';
    if (fecs.n_elem > N) Rcout << "... ";
    Rcout << fecs(fecs.n_elem-1) << ')' << endl;
    
    
    Rcout << "Population numbers:" << endl;
    
    Rcout << "  - aphids[t]:   (";
    for (unsigned i = 0; i < (N-1); i++) Rcout << X_t(i) << ' ';
    if (X_t.n_elem > N) Rcout << "... ";
    Rcout << X_t(X_t.n_elem-1) << ')' << endl;
    
    Rcout << "  - aphids[t+1]: (";
    for (unsigned i = 0; i < (N-1); i++) Rcout << X_t1(i) << ' ';
    if (X_t1.n_elem > N) Rcout << "... ";
    Rcout << X_t1(X_t1.n_elem-1) << ')' << endl;
    
    N = std::min(static_cast<arma::uword>(6), Y_t.n_elem);
    Rcout << "  - wasps[t]:    (";
    for (unsigned i = 0; i < (N-1); i++) Rcout << Y_t(i) << ' ';
    if (Y_t.n_elem > N) Rcout << "... ";
    Rcout << Y_t(Y_t.n_elem-1) << ')' << endl;
    
    Rcout << "  - wasps[t+1]:  (";
    for (unsigned i = 0; i < (N-1); i++) Rcout << Y_t1(i) << ' ';
    if (Y_t.n_elem > N) Rcout << "... ";
    Rcout << Y_t1(Y_t1.n_elem-1) << ')' << endl;
    
    return;
}

// Defining adding process error
void AphidWasp::process_error(double z, double Y_m) {
    
    if (demog_mult == 0 && sigma_x == 0 && sigma_y == 0) return;
    
    uint living_aphids = X_t.n_elem + mum_days(0);
    uint total_stages = X_t.n_elem + Y_t.n_elem;
    
    arma::mat Se(total_stages, total_stages, arma::fill::zeros);
    
    // Aphid (both parasitized and not) process error
    Se(arma::span(0, living_aphids-1), arma::span(0, living_aphids-1)) =
        (sigma_x*sigma_x + demog_mult * std::min(0.5, 1 / std::abs(1 + z))) *
        (rho * arma::mat(living_aphids,living_aphids,arma::fill::ones) + 
        (1-rho) * arma::mat(living_aphids,living_aphids,arma::fill::eye));
    
    // Mummy process error needs no action because it's already zero

    // Adult parasitoid process error
    Se(Se.n_rows-1, Se.n_cols-1) = sigma_y*sigma_y + 
        demog_mult * std::min(0.5, 1 / std::abs(1 + Y_m));
    
    // chol doesn't work with zeros on diagonal
    arma::uvec non_zero = arma::find(Se.diag() > 0);
    
    // Cholesky decomposition of Se so output has correct variance-covariance matrix
    //   "a vector of independent normal random variables,
    //   when multiplied by the transpose of the Cholesky deposition of [Se] will
    //   have covariance matrix equal to [Se]."
    arma::mat chol_decomp = arma::chol(Se(non_zero,non_zero)).t();
    
    // Random numbers from distribution N(0,1)
    arma::vec E(non_zero.n_elem);
    for (uint i = 0; i < non_zero.n_elem; i++) E(i) = rnorm_distr(rnd_engine);
    
    // Making each element of E have correct variance-covariance matrix
    E = chol_decomp * E;
    
    arma::uvec nz_aphid = non_zero(arma::find(non_zero < X_t1.n_rows));
    arma::uvec nz_wasp = non_zero(arma::find(non_zero >= X_t1.n_rows)) - X_t1.n_rows;
    
    // Plugging in errors into the X[t+1] and Y[t+1] vectors
    X_t1.rows(nz_aphid) = X_t1.rows(nz_aphid) % exp(E.head(nz_aphid.n_elem));
    Y_t1.rows(nz_wasp) = Y_t1.rows(nz_wasp) % exp(E.tail(nz_wasp.n_elem));
    
    // Because we used normal distributions to approximate demographic and environmental 
    // stochasticity, it is possible for aphids and parasitoids to 
    // "spontaneously appear" when the estimate of e(t) is large. To disallow this 
    // possibility, the number of aphids and parasitized aphids in a given age class 
    // on day t was not allowed to exceed the number in the preceding age class on 
    // day t – 1.
    
    for (uint i = 1; i < X_t1.n_elem; i++) {
        if (X_t1(i) > X_t(i-1)) X_t1(i) = X_t(i-1);
    }
    // Not going to the end for parasitoids bc you can have more adults than mummies
    // bc adults stay in that stage for >1 days
    for (uint i = 1; i < (Y_t1.n_elem-1); i++) {
        if (Y_t1(i) > Y_t(i-1)) Y_t1(i) = Y_t(i-1);
    }
    
    return;
}


// Calculate X[t+1]
void AphidWasp::iterate_X(double S) {
    X_t = X_t1;
    X_t1 = (pred_rate * S * A) % (leslie * X_t);
    return;
}

// Calculate Y(t+1)
void AphidWasp::iterate_Y(double S_y) {
    
    Y_t = Y_t1;
    
    const arma::mat& L(leslie);
    double m_1 = mum_days(0);
    
    arma::mat tmp;
    // Y_1(t+1)                             // initial parasitized aphids
    tmp = S_y * (1 - A).t() * (L * X_t);
    // below, arma::as_scalar is to let this know it's a double; it's already 1x1
    Y_t1(0) = arma::as_scalar(tmp);
    arma::vec s_i = arma::diagvec(L, -1);
    // Y_i(t+1) for (i = 1, ..., m_1)       // parasitized but alive aphids
    Y_t1(arma::span(1, m_1)) = (s_i.head(m_1) * S_y) % Y_t.head(m_1);
    // Y_i(t+1) for (i = m_1+1, ..., m-1)   // mummies
    Y_t1(arma::span(m_1+1, Y_t1.n_elem-2)) = pred_rate * 
        Y_t(arma::span(m_1, Y_t.n_elem-3));
    // Y_m(t+1)                             // adult wasps
    Y_t1.tail(1) = s_y * Y_t.tail(1) + sex_ratio * Y_t(Y_t.n_elem-2);
    
    return;
}










// ======================================================================================
// ======================================================================================
// ======================================================================================

//          One patch, multiple aphid lines

// ======================================================================================
// ======================================================================================
// ======================================================================================





// class OnePatch
    
// Density dependence for aphids (note that equation is different in paper)
double OnePatch::S(double K) {
    return 1 / (1 + K * z);
};
    
// Density dependence for wasps (note that equation is different in paper)
double OnePatch::S_y(double K_y) {
    return 1 / (1 + K_y * z);
};
    


void OnePatch::show() {
    Rcout << "< Info for one patch w " << pops.size() << " aphid-wasp combos >" << endl;
    // Rcout << "harvest period: " << harvest_period << endl;
    // Rcout << "harvest offset: " << harvest_offset << endl;
    Rcout << "Current numbers:" << endl;
    Rcout << "  - non-parasitized aphids: " << x << endl;
    Rcout << "  - living aphids:          " << z << endl;
    Rcout << "  - adult wasps:            " << Y_m << endl;
}


// Boolean for whether to harvest at time t
bool OnePatch::do_harvest(uint t) {
    if (t < harvest_offset || harvest_period == 0 || t == 0) return false;
    uint a = t - harvest_offset;
    uint b = a % harvest_period;
    return b == 0;
};


// iterate patch, doing everything but dispersal
void OnePatch::iterate_patch(uint t) {
    
    // Reset, then recalculate patch totals
    z = 0;
    x = 0;
    Y_m = 0;
    for (AphidWasp& aw : pops) {
        double ta = aw.total_aphids();
        double tlp = aw.total_living_paras();
        z += ta + tlp;
        x += ta;
        Y_m += aw.total_adult_wasps();
    }
    
    // Update each AphidWasp object in pops
    for (AphidWasp& aw : pops) {
        // Equivalent to S(z(t)) and S_y(z(t))
        double St = S(aw.K);
        double S_yt = S_y(aw.K_y);
        // Update attack probabilities
        aw.iterate_A(Y_m, x);
        // Update aphid and wasp densities
        aw.iterate_X(St);
        aw.iterate_Y(S_yt);
        aw.process_error(z, Y_m);
        bool harvest_bool = do_harvest(t);
        if (harvest_bool) aw.harvest();
    }
    return;
}


// Update output summary for one time, one line, one patch
void OnePatch::update_summary(SimSummary& output, uint t_, uint p_) {
    
    for (uint l_ = 0; l_ < pops.size(); l_++) {
        AphidWasp& aw(pops[l_]);
        output.fill(
            t_, l_, p_,
            arma::sum(aw.X_t1),  // aphids
            arma::sum(aw.Y_t1(arma::span(0,aw.mum_days(0)-1))),  // parasit
            arma::sum(aw.Y_t1(arma::span(aw.mum_days(0),aw.Y_t1.n_elem-2))),  // mummies
            arma::as_scalar(aw.Y_t1.tail(1))  // wasps
        );
    }
    return;
}

// Reset and set new seed if you want to do another simulation set
void OnePatch::reset_patch(uint rng_seed) {
    sitmo::prng_engine eng(rng_seed);
    z = 0;
    x = 0;
    Y_m = 0;
    for (AphidWasp& aw : pops) {
        aw.A = arma::zeros<arma::vec>(aw.rel_attack.n_elem);
        aw.Y_t = aw.Y_0;
        aw.Y_t1 = aw.Y_0;
        aw.X_t = aw.X_0;
        aw.X_t1 = aw.X_0;
        aw.set_seed(eng());
    }
}


// Return reference to an AphidWasp object of a given name
AphidWasp& OnePatch::find_line(string aphid_name_) {
    for (AphidWasp& aw : pops) {
        if (aw.aphid_name == aphid_name_) return aw;
    }
    stop("Aphid name not found");
}







// ======================================================================================
// ======================================================================================
// ======================================================================================

//          All Patches

// ======================================================================================
// ======================================================================================
// ======================================================================================


// Dispersal among patches
void SimPatches::dispersal() {
    // Do this individually for each aphid-wasp combo
    for (const string& an : aphid_names) {
        // Calculate the mean number of aphids (for each stage) and 
        // wasps (just adults) across patches
        arma::vec mean_aphids = arma::zeros<arma::vec>(99); // 99 is more than enough
        double mean_wasps = 0;
        double n_patches = 0;
        for (OnePatch& p : patches) {
            AphidWasp& aw = p.find_line(an);
            mean_aphids.head(aw.n_aphid_stages) += aw.X_t1;
            mean_wasps += arma::as_scalar(aw.Y_t1.tail(1));
            n_patches++;
        }
        mean_aphids /= n_patches;
        mean_wasps /= n_patches;
        
        // Now go back through and actually carry out the dispersal
        for (OnePatch& p : patches) {
            AphidWasp& aw = p.find_line(an);
            /* Aphid dispersal */
            // Indices (hence the -1) for stages at which dispersal occurs
            arma::uvec disp_stages = arma::regspace<arma::uvec>(aw.disp_start-1, 
                                                                aw.n_aphid_stages-1);
            double immigration;
            for (arma::uword i : disp_stages) {
                immigration = aw.disp_aphid * mean_aphids(i);
                aw.X_t1(i) *= (1 - aw.disp_aphid); // emigration
                aw.X_t1(i) += immigration;
            }
            /* Wasp dispersal */
            immigration = aw.disp_wasp * mean_wasps;
            aw.Y_t1.tail(1) *= (1 - aw.disp_wasp); // emigration
            aw.Y_t1.tail(1) += immigration;
        }
    }
    return;
}


// reset object and simulate a set number of time periods
SimSummary SimPatches::simulate(uint max_t) {
    
    uint rng_seed = static_cast<uint>(R::runif(0, MAX_INT));
    
    sitmo::prng_engine eng(rng_seed);
    // Make sure everything is back to normal and new, unique seeds are set
    for (OnePatch& p : patches) p.reset_patch(eng());
    // Reset time to zero
    t = 0;
    
    SimSummary output(max_t, aphid_names.size(), patches.size());
    
    while (t < max_t) {
        for (OnePatch& p : patches) p.iterate_patch(t);
        dispersal();
        for (uint p_ = 0; p_ < patches.size(); p_++) {
            OnePatch& p(patches[p_]);
            p.update_summary(output, t, p_);
        }
        t++;
    }
    
    return output;
}


