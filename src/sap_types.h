# ifndef _SAP_TYPES_H
# define _SAP_TYPES_H

#include <RcppArmadillo.h> // arma namespace
#include <sitmo.h>         // parallel rng
#include <vector>          // vector class
#include <cmath>           // log, exp
#include <random>          // normal distribution
#include <cstdint>         // integer types
#include <algorithm>       // find
#include "math.h"          // leslie_matrix and leslie_sad


using namespace Rcpp;
using namespace std;


typedef uint_fast8_t uint8;
typedef uint_fast32_t uint;
typedef int_fast32_t sint;
typedef uint_fast64_t uint64;
typedef int_fast64_t sint64;





// Summarize a simulation's populations for all lines across all times and patches
// Rows are time, columns are lines, slices (3rd dimension) are patches
//' @export SimSummary
struct SimSummary {

    arma::cube aphids;          // density of unparasitized aphids
    arma::cube parasit;         // density of parasitized, but alive, aphids
    arma::cube mummies;         // density of mummies
    arma::cube wasps;           // density of adult wasps
    
    // ------------------
    // Constructors
    // ------------------
    
    SimSummary() : aphids(), parasit(), mummies(), wasps() {};
    SimSummary(uint max_t, uint n_lines, uint n_patches) 
        : aphids(max_t, n_lines, n_patches), 
          parasit(max_t, n_lines, n_patches),
          mummies(max_t, n_lines, n_patches), 
          wasps(max_t, n_lines, n_patches) {};
    
    
    // ------------------
    // Overloaded fill function:
    // ------------------
    
    // For one time, one line, one patch
    void fill(uint t_, uint l_, uint p_, 
              double aphids_, double parasit_, double mummies_, double wasps_) {
        aphids(t_, l_, p_) = aphids_;
        parasit(t_, l_, p_) = parasit_;
        mummies(t_, l_, p_) = mummies_;
        wasps(t_, l_, p_) = wasps_;
    }
    
    // one time, one patch, ALL lines
    void fill(uint t_, uint p_, 
              arma::rowvec aphids_, arma::rowvec parasit_, 
              arma::rowvec mummies_, arma::rowvec wasps_) {
        aphids.slice(p_).row(t_) = aphids_;
        parasit.slice(p_).row(t_) = parasit_;
        mummies.slice(p_).row(t_) = mummies_;
        wasps.slice(p_).row(t_) = wasps_;
    }
    
    // one time, ALL patches, ALL lines
    // For the input matrices, rows are patches, columns are lines
    void fill(uint t_,
              arma::mat aphids_, arma::mat parasit_,
              arma::mat mummies_, arma::mat wasps_) {
        for (uint i = 0; i < aphids.n_slices; i++) {
            aphids.slice(i).row(t_) = aphids_.row(i);
            parasit.slice(i).row(t_) = parasit_.row(i);
            mummies.slice(i).row(t_) = mummies_.row(i);
            wasps.slice(i).row(t_) = wasps_.row(i);
        }
    }
    
};




// ======================================================================================
// ======================================================================================
// ======================================================================================

//          One aphid line and associated wasps

// ======================================================================================
// ======================================================================================
// ======================================================================================

// Aphid population
struct AphidPop {
    
    const arma::mat leslie;       // Leslie matrix with survival and reproduction
    const arma::vec X_0;          // initial aphid abundances by stage
    const double K;               // aphid density dependence
    const uint n_aphid_stages;    // number of aphid stages (i.e., days)
    
    // Changing through time
    arma::vec X_t;                // Aphid density at time t
    arma::vec X_t1;               // Aphid density at time t+1
    
    // Constructors
    AphidPop(List par_list) 
        : leslie(leslie_matrix(as<arma::uvec>(par_list["instar_days"]),
                 as<double>(par_list["surv_juv"]),
                 as<arma::vec>(par_list["surv_adult"]),
                 as<arma::vec>(par_list["repro"]))),
         X_0(as<double>(par_list["aphid_density_0"]) * leslie_sad(leslie)), 
         K(as<double>(par_list["K"])), 
         n_aphid_stages(arma::sum(as<arma::uvec>(par_list["instar_days"]))),
         X_t(X_0), X_t1(X_0) {};
    
    // Total (non-parasitized) aphids
    double total_aphids() {
        return arma::sum(X_t1);
    }
    
};


// Wasp population
struct WaspPop {
    
    const arma::vec Y_0;          // initial wasp abundances by stage
    const double sex_ratio;       // proportion of female wasps
    const double K_y;             // parasitized aphid density dependence
    const double s_y;             // parasitoid adult daily survival
    const arma::uvec mum_days;    // number of days per mummy stage (aphid alive & dead)
    const uint n_wasp_stages;     // number of wasp stages (i.e., days)
    
    // Changing through time
    arma::vec Y_t;                // Wasp density at time t
    arma::vec Y_t1;               // Wasp density at time t+1
    
    // Constructor
    WaspPop(List par_list) 
        : Y_0(arma::join_cols(
                arma::zeros<arma::vec>(arma::sum(as<arma::uvec>(par_list["mum_days"]))),
                as<double>(par_list["wasp_density_0"]) * arma::ones<arma::vec>(1))),
          sex_ratio(as<double>(par_list["sex_ratio"])), 
          K_y(as<double>(par_list["K_y"])), 
          s_y(as<double>(par_list["s_y"])), 
          mum_days(as<arma::uvec>(par_list["mum_days"])),
          n_wasp_stages(arma::sum(mum_days) + 1),
          Y_t(Y_0), Y_t1(Y_0) {};

    // Total living, but parasitized aphids
    double total_living_paras() {
        return arma::sum(Y_t1(arma::span(0, mum_days(0)-1)));
    }
    
    // Total adult wasps
    double total_adult_wasps() {
        return Y_t1(Y_t1.n_elem-1);
    }

};






// Wasp attack
struct WaspAttack {
    
    const arma::vec rel_attack;   // relative wasp attack rates by aphid stage
    const double a;               // overall parasitoid attack rate
    const double k;               // aggregation parameter of the nbinom distribution
    const double h;               // parasitoid attack rate handling time
    const arma::vec attack_surv;  // survival rates of singly & multiply attacked aphids
    
    // Changing through time
    arma::vec A;                  // attack probabilities at time t
    
    // Constructor
    WaspAttack(List par_list) 
        : rel_attack(as<arma::vec>(par_list["rel_attack"])),
          a(as<double>(par_list["a"])), 
          k(as<double>(par_list["k"])), 
          h(as<double>(par_list["h"])), 
          attack_surv(as<arma::vec>(par_list["attack_surv"])),
          A(arma::zeros<arma::vec>(rel_attack.n_elem)) {};
    
    // 
    // Update attack probabilities
    // Equation 6 from Meisner et al. (2014)
    // Note: rel_attack is equivalent to p_i
    // 
    void iterate_A(double Y_m, double x) {

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
    
};


// Process error
struct ProcessError {

    const double sigma_x;         // environmental standard deviation for aphids
    const double sigma_y;         // environmental standard deviation for parasitoids
    const double rho;             // environmental correlation among instars
    const double demog_mult;      // multiplier for demographic stochasticity
    
    // Constructor
    ProcessError(List par_list) 
        : sigma_x(as<double>(par_list["sigma_x"])), 
          sigma_y(as<double>(par_list["sigma_y"])), 
          rho(as<double>(par_list["rho"])), 
          demog_mult(as<double>(par_list["demog_mult"])) {};

};



// Environment: harvest, dispersal, and predation
struct Environ {

    const double harvest_surv;    // survival rate for living aphids during a harvest
    const double disp_aphid;      // dispersal rate for aphids
    const double disp_wasp;       // dispersal rate for wasps
    const uint disp_start;        // stage in which dispersal starts in aphids
    const double pred_rate;       // predation on aphids and mummies
    
    // Constructor
    Environ(List par_list) 
        : harvest_surv(as<double>(par_list["harvest_surv"])),
          disp_aphid(as<double>(par_list["disp_aphid"])),
          disp_wasp(as<double>(par_list["disp_wasp"])),
          disp_start(as<uint>(par_list["disp_start"])),
          pred_rate(as<double>(par_list["pred_rate"])) {};
    
};






// Info about one aphid line and wasps that parasitize them

//' @export AphidWasp
class AphidWasp: public AphidPop, public WaspPop, public WaspAttack, 
                 public ProcessError, public Environ  {
public:

    // --------
    // Members:
    // --------

    const string aphid_name;    // unique identifying name for this aphid line

    // --------
    // Constructor:
    // --------
    
    AphidWasp(string aphid_name_, List par_list, uint seed)
        : AphidPop::AphidPop(par_list), 
          WaspPop::WaspPop(par_list), 
          WaspAttack::WaspAttack(par_list), 
          ProcessError::ProcessError(par_list), 
          Environ::Environ(par_list),
          aphid_name(aphid_name_),
          rnorm_distr(0,1),
          rnd_engine(seed) {};
    
    // ----------------
    // Methods (all but set_seed definitions below)
    // ----------------
    
    // Show an object in cmd
    void show() const;
    
    // Add process error
    void process_error(double z, double Y_m);
    
    void iterate_X(double S);
    
    void iterate_Y(double S_y);
    
    void harvest();
    
    void set_seed(uint seed) {
        rnd_engine.seed(seed);
    }

    // ----------------
    // Private members (for rng):
    // ----------------
private:
    std::normal_distribution<double> rnorm_distr;
    sitmo::prng_engine rnd_engine;
};



void AphidWasp::harvest() {
    // Kill non-parasitized aphids
    X_t1 *= harvest_surv;
    // Kill parasitized (but still living) aphids at the same rate
    Y_t1.head(mum_days(0)) *= harvest_surv;
    // Kill all mummies
    Y_t1(arma::span(mum_days(0), Y_t1.n_elem-2)).fill(0);
}






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
    
    // Below is not needed since they're already zeros
    // // Mummy process error, turning back to zero
    // Se(arma::span(living_aphids,Se.n_rows-2), 
    //    arma::span(living_aphids,Se.n_cols-2)).fill(0);
    // Se(arma::span(0,living_aphids-1), 
    //    arma::span(living_aphids,Se.n_cols-2)).fill(0);
    // Se(arma::span(living_aphids,Se.n_rows-2), 
    //    arma::span(0,living_aphids-1)).fill(0);
    
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





class OnePatch {
public:
    
    // -------
    // Members:
    // -------
    const uint harvest_period;        // time points between harvests
    const uint harvest_offset;        // time at which to begin harvests
    double z;                         // Sum of all living aphids at time t
    double x;                         // Sum of non-parasitized aphids at time t
    double Y_m;                       // Total number of adult wasps
    vector<AphidWasp> pops;           // Vector of aphid-wasp combos
    
    // -------
    // Constructor:
    // -------
    
    OnePatch(vector<string> aphid_names, vector<List> par_L, vector<uint> seeds,
             const uint& harvest_period_, const uint& harvest_offset_)
        : harvest_period(harvest_period_), 
          harvest_offset(harvest_offset_),
          z(0.0), x(0.0), Y_m(0.0) {
        
        if (aphid_names.size() != par_L.size() || par_L.size() != seeds.size()) {
            stop("aphid_names, par_L, and seeds must have the same length.");
        }
        for (uint i = 0; i < aphid_names.size(); i++) {
            AphidWasp aw(aphid_names[i], par_L[i], seeds[i]);
            pops.push_back(aw);
        }
    };
    
    
    // Methods:
    
    // Density dependence for aphids (note that equation is different in paper)
    double S(double K) {
        return 1 / (1 + K * z);
    };
    
    // Density dependence for wasps (note that equation is different in paper)
    double S_y(double K_y) {
        return 1 / (1 + K_y * z);
    };
    
    // Boolean for whether to harvest at time t
    bool do_harvest(uint t);
    
    // Return reference to an AphidWasp object of a given name
    AphidWasp& find_line(string aphid_name_);
    
    // iterate patch, doing everything but dispersal
    void iterate_patch(uint t);
    
    // retrieve info from a patch (to be done after dispersal on SimPatches)
    // and update a SimSummary object
    void update_summary(SimSummary& output, uint t_, uint p_);

    // Reset and set new seed if you want to do another simulation set
    void reset_patch(uint rng_seed);
    
};


// Boolean for whether to harvest at time t
bool OnePatch::do_harvest(uint t) {
    if (t < harvest_offset) return false;
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
        // Equivalent to S(z(t)) and S_y(z(t)) [no idea why equations are different]
        double St = S(aw.K);
        double S_yt = S_y(aw.K_y);
        // Update attack probabilities
        aw.iterate_A(Y_m, x);
        // Update aphid and wasp densities
        aw.iterate_X(St);
        aw.iterate_Y(S_yt);
        aw.process_error(z, Y_m);
        bool do_harvest_ = do_harvest(t);
        if (do_harvest_) aw.harvest();
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









//' @export SimPatches
class SimPatches {
public:
    
    // -------
    // Members:
    // -------
    
    const vector<string> aphid_names;   // vector of aphid names (same for all patches)
    vector<OnePatch> patches;           // vector of patch info
    uint t;                             // current time
    
    // -------
    // Constructor:
    // -------
    
    // aphid_names is not nested bc all patches should have the same names 
    // I want them to be allowed to have different parameters, though, to
    // simulate environmental effects
    SimPatches(const vector<string>& aphid_names_,
               const vector<vector<List>>& par_L, 
               const uint& seed,
               const vector<uint>& harvest_periods, 
               const vector<uint>& harvest_offsets) 
        : aphid_names(aphid_names_), t(0) {
        
        uint n_patches = harvest_periods.size();
        
        sitmo::prng_engine eng(seed);
        
        if (n_patches != harvest_offsets.size()) {
            stop("harvest_periods and harvest_offsets should have the same length.");
        }
        for (uint i = 0; i < n_patches; i++) {
            vector<uint> seeds(par_L.size());
            for (uint j = 0; j < par_L[i].size(); j++) seeds[j] = eng();
            OnePatch op(aphid_names, par_L[i], seeds,
                        harvest_periods[i], harvest_offsets[i]);
            patches.push_back(op);
        }
    };
    
    void dispersal();
    
    SimSummary simulate(uint max_t, uint rng_seed);
    
};




void SimPatches::dispersal() {
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
            arma::uvec disp_stages = arma::regspace<arma::uvec>(aw.disp_start-1, 
                                                                aw.n_aphid_stages-1);
            double immigration;
            for (arma::uword i : disp_stages) {
                immigration = aw.disp_aphid * mean_aphids(i);
                aw.X_t1(i) *= (1 - aw.disp_aphid); // emigration
                aw.X_t1(i) += immigration;
            }
            immigration = aw.disp_wasp * mean_wasps;
            aw.Y_t1.tail(1) *= (1 - aw.disp_wasp); // emigration
            aw.Y_t1.tail(1) += immigration;
        }
    }
    return;
}


SimSummary SimPatches::simulate(uint max_t, uint rng_seed) {
    
    sitmo::prng_engine eng(rng_seed);
    // Make sure everything is back to normal and new, unique seeds are set
    for (OnePatch& p : patches) p.reset_patch(eng());
    
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




#endif