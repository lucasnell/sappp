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





// Aphid population
struct aphid_pop {
    
    const arma::mat leslie;       // Leslie matrix with survival and reproduction
    const arma::vec X_0;          // initial aphid abundances by stage
    const double K;               // aphid density dependence
    const uint n_stages;          // number of aphid stages (i.e., days)
    
    // Changing through time
    arma::vec X_t;                // Aphid density at time t
    arma::vec X_t1;               // Aphid density at time t+1
    
    // Constructors
    aphid_pop(List par_list) 
        : leslie(leslie_matrix(as<arma::uvec>(par_list["instar_days"]),
                 as<double>(par_list["surv_juv"]),
                 as<arma::vec>(par_list["surv_adult"]),
                 as<arma::vec>(par_list["repro"]))),
         X_0(as<double>(par_list["aphid_density_0"]) * leslie_sad(leslie)), 
         K(as<double>(par_list["K"])), 
         n_stages(arma::sum(as<arma::uvec>(par_list["instar_days"]))),
         X_t(X_0), X_t1(X_0) {};
    
};


// Wasp population
struct wasp_pop {
    
    // I'm only putting mum_days here so it gets initialized first
    const arma::vec Y_0;          // initial wasp abundances by stage
    const double sex_ratio;       // proportion of female wasps
    const double K_y;             // parasitized aphid density dependence
    const double s_y;             // parasitoid adult daily survival
    const arma::uvec mum_days;    // number of days per mummy stage (aphid alive & dead)
    const uint n_stages;          // number of wasp stages (i.e., days)
    
    // Changing through time
    arma::vec Y_t;                // Wasp density at time t
    arma::vec Y_t1;               // Wasp density at time t+1
    
    // Constructor
    wasp_pop(List par_list) 
        : Y_0(arma::join_cols(
                arma::zeros<arma::vec>(arma::sum(as<arma::uvec>(par_list["mum_days"]))),
                as<double>(par_list["wasp_density_0"]) * arma::ones<arma::vec>(1))),
          sex_ratio(as<double>(par_list["sex_ratio"])), 
          K_y(as<double>(par_list["K_y"])), 
          s_y(as<double>(par_list["s_y"])), 
          mum_days(as<arma::uvec>(par_list["mum_days"])),
          n_stages(arma::sum(mum_days) + 1) {};

};

// Wasp attack
struct wasp_attack {
    
    const arma::vec rel_attack;   // relative wasp attack rates by aphid stage
    const double a;               // overall parasitoid attack rate
    const double k;               // aggregation parameter of the nbinom distribution
    const double h;               // parasitoid attack rate handling time
    const arma::vec attack_surv;  // survival rates of singly & multiply attacked aphids
    
    // Changing through time
    arma::vec A;                  // attack probabilities at time t
    
    // Constructor
    wasp_attack(List par_list) 
        : rel_attack(as<arma::vec>(par_list["rel_attack"])),
          a(as<double>(par_list["a"])), 
          k(as<double>(par_list["k"])), 
          h(as<double>(par_list["h"])), 
          attack_surv(as<arma::vec>(par_list["attack_surv"])),
          A(arma::zeros<arma::vec>(rel_attack.n_elem)) {};
    
};


// Process error
struct process_error {
    const double sigma_x;         // environmental standard deviation for aphids
    const double sigma_y;         // environmental standard deviation for parasitoids
    const double rho;             // environmental correlation among instars
    const double demog_mult;      // multiplier for demographic stochasticity
    
    // Constructor
    process_error(List par_list) 
        : sigma_x(as<double>(par_list["sigma_x"])), 
          sigma_y(as<double>(par_list["sigma_y"])), 
          rho(as<double>(par_list["rho"])), 
          demog_mult(as<double>(par_list["demog_mult"])) {};

};



// Environment: harvest, dispersal, and predation
struct environ {
    const double harvest_surv;    // survival rate for living aphids during a harvest
    const double disp_aphid;      // dispersal rate for aphids
    const double disp_wasp;       // dispersal rate for wasps
    const uint disp_start;        // stage in which dispersal starts in aphids
    const double pred_rate;       // predation on aphids and mummies
    
    // Constructor
    environ(List par_list) 
        : harvest_surv(as<double>(par_list["harvest_surv"])),
          disp_aphid(as<double>(par_list["disp_aphid"])),
          disp_wasp(as<double>(par_list["disp_wasp"])),
          disp_start(as<uint>(par_list["disp_start"])),
          pred_rate(as<double>(par_list["pred_rate"])) {};
    
};






// Info about one aphid line and wasps that parasitize them

//' @export aphid_wasp
class aphid_wasp {
public:

    // --------
    // Members:
    // --------

    const string aphid_name;    // unique identifying name for this aphid line
    
    // Aphid population
    aphid_pop aphids;

    // Wasp population
    wasp_pop wasps;

    // Wasp attack
    wasp_attack attacks;

    // Process error
    process_error errors;

    // Environment: harvest, dispersal, and predation
    environ envir;


    // --------
    // Constructors:
    // --------
    
    aphid_wasp(List par_list) 
        : aphid_name(as<string>(par_list["aphid_name"])), 
          aphids(par_list), wasps(par_list), attacks(par_list), errors(par_list),
          envir(par_list) {};

    void show() const {
        arma::rowvec survs = arma::diagvec(aphids.leslie, -1).t();
        arma::rowvec fecs = aphids.leslie(0, arma::span(1, aphids.leslie.n_cols - 1));
        
        arma::uword N = std::min(static_cast<arma::uword>(6), aphids.leslie.n_rows);
        
        Rcout.precision(4);
        Rcout << std::fixed;
        
        Rcout << "Constant info for '" << aphid_name << "' line" << endl;
        Rcout << "  * Resistances: (";
        Rcout << attacks.attack_surv(0) << ' ' << attacks.attack_surv(1) << ')' << endl;
        
        Rcout << "  * Survivals:   (";
        for (unsigned i = 0; i < N; i++) Rcout << survs(i) << ' ';
        if (aphids.leslie.n_rows > N) Rcout << "...";
        Rcout << ')' << endl;
        
        Rcout << "  * Fecundities: (";
        for (unsigned i = 0; i < N; i++) Rcout << fecs(i) << ' ';
        if (aphids.leslie.n_rows > N) Rcout << "...";
        Rcout << ')' << endl;
                
        return;
    }

};









//' @export patch
class patch {
public:
    
    // -------
    // Members:
    // -------
    uint harvest_period;       // time points between harvests
    uint harvest_offset;       // time at which to begin harvests
    double z;                  // Sum of all living aphids at time t
    double x;                  // Sum of non-parasitized aphids at time t
    double Y_m;                // Total number of adult wasps
    
    // -------
    // Constructor:
    // -------
    
    patch(const uint& harvest_period_, const uint& harvest_offset_)
        : harvest_period(harvest_period_), harvest_offset(harvest_offset_),
          z(0.0), x(0.0), Y_m(0.0) {};
    
    
    // Methods:
    
    // Density dependence for aphids (note that equation is different in paper)
    double S(double K) {
        return 1 / (1 + K * z);
    };
    
    // Density dependence for wasps (note that equation is different in paper)
    double S_y(double K_y) {
        return 1 / (1 + K_y * z);
    };
    
    
    
};


#endif