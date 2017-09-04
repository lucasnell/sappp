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
class SimSummary {
public:

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
    
    
    /* Overloaded fill function: */
    
    // For one time, one line, one patch
    void fill(uint t_, uint l_, uint p_, 
              double aphids_, double parasit_, double mummies_, double wasps_);
    // one time, one patch, ALL lines
    void fill(uint t_, uint p_, 
              arma::rowvec aphids_, arma::rowvec parasit_, 
              arma::rowvec mummies_, arma::rowvec wasps_);
    // one time, ALL patches, ALL lines
    // For the input matrices, rows are patches, columns are lines
    void fill(uint t_,
              arma::mat aphids_, arma::mat parasit_,
              arma::mat mummies_, arma::mat wasps_);
    
    void show();
    
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
    double total_aphids();
    
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
    double total_living_paras();
    
    // Total adult wasps
    double total_adult_wasps();

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
    void iterate_A(double Y_m, double x);
    
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

class AphidWasp: public AphidPop, public WaspPop, public WaspAttack, 
                 public ProcessError, public Environ  {
public:

    // --------
    // Members:
    // --------

    const string aphid_name;    // unique identifying name for this aphid-wasp combo

    // --------
    // Constructor:
    // --------
    
    AphidWasp(List par_list)
        : AphidPop::AphidPop(par_list), 
          WaspPop::WaspPop(par_list), 
          WaspAttack::WaspAttack(par_list), 
          ProcessError::ProcessError(par_list), 
          Environ::Environ(par_list),
          aphid_name(as<string>(par_list["aphid_name"])),
          rnorm_distr(0,1),
          rnd_engine() {};
    
    // ----------------
    // Methods
    // ----------------
    
    // Show an object in cmd
    void show() const;
    
    // Add process error
    void process_error(double z, double Y_m);
    
    void iterate_X(double S);
    
    void iterate_Y(double S_y);
    
    void harvest();
    
    void set_seed(uint seed);

    // ----------------
    // Private members (for rng):
    // ----------------
private:
    std::normal_distribution<double> rnorm_distr;
    sitmo::prng_engine rnd_engine;
};






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
    const uint harvest_period;          // time points between harvests
    const uint harvest_offset;          // time at which to begin harvests
    double z;                           // Sum of all living aphids at time t
    double x;                           // Sum of non-parasitized aphids at time t
    double Y_m;                         // Total number of adult wasps
    vector<AphidWasp> pops;             // Vector of aphid-wasp combos
    
    // -------
    // Constructor:
    // -------
    
    OnePatch(vector<List> par_L, 
             const uint& harvest_period_,
             const uint& harvest_offset_)
         : harvest_period(harvest_period_),
           harvest_offset(harvest_offset_),
          z(0.0), x(0.0), Y_m(0.0) {
               
        for (uint i = 0; i < par_L.size(); i++) {
            AphidWasp aw(par_L[i]);
            pops.push_back(aw);
        }
    };
    
    
    // Methods:
    
    void show();
    
    // Density dependence for aphids (note that equation is different in paper)
    double S(double K);
    
    // Density dependence for wasps (note that equation is different in paper)
    double S_y(double K_y);
    
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










// ======================================================================================
// ======================================================================================
// ======================================================================================

//          All Patches

// ======================================================================================
// ======================================================================================
// ======================================================================================





//' @export SimPatches
class SimPatches {
public:
    
    // -------
    // Members:
    // -------
    
    vector<string> aphid_names;         // vector of aphid names (same for all patches)
    vector<OnePatch> patches;           // vector of patch info
    uint t;                             // current time
    
    // -------
    // Constructor:
    // -------
    
    // aphid_names is not nested bc all patches should have the same names 
    // I want them to be allowed to have different parameters, though, to
    // simulate environmental effects
    SimPatches(vector<vector<List>> par_L,
               vector<uint> harvest_periods,
               vector<uint> harvest_offsets)
        : aphid_names(0), t(0) {
        
        
        uint n_patches = harvest_periods.size();
        
        if (n_patches != harvest_offsets.size()) {
            stop("harvest_periods and harvest_offsets should have the same length.");
        }
        
        // Going through first patch's info to get line names
        // They should be the same among patches
        uint n_pops = par_L[0].size();
        for (uint i = 0; i < n_pops; i++) {
            aphid_names.push_back(as<string>(par_L[0][i]["aphid_name"]));
        }
        
        for (uint i = 0; i < n_patches; i++) {
            OnePatch op(par_L[i], harvest_periods[i], harvest_offsets[i]);
            // Checking that names are the same
            for (uint j = 0; j < op.pops.size(); j++) {
                if (op.pops[j].aphid_name != aphid_names[j]) {
                    stop("All patches' aphid names should be identical.");
                }
            }
            patches.push_back(op);
        }
    };
    
    void dispersal();
    
    SimSummary simulate(uint max_t, uint rng_seed);
    
};




#endif