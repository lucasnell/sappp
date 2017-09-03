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
class AphidWasp: public AphidPop, WaspPop, WaspAttack, ProcessError, Environ  {
public:

    // --------
    // Members:
    // --------

    const string aphid_name;    // unique identifying name for this aphid line
    
    // // Aphid population
    // AphidPop aphids;
    // 
    // // Wasp population
    // WaspPop wasps;
    // 
    // // Wasp attack
    // WaspAttack attacks;
    // 
    // // Process error
    // ProcessError errors;
    // 
    // // Environment: harvest, dispersal, and predation
    // Environ envir;


    // --------
    // Constructors:
    // --------
    
    AphidWasp(string aphid_name_, List par_list)
        : AphidPop::AphidPop(par_list), 
          WaspPop::WaspPop(par_list), 
          WaspAttack::WaspAttack(par_list), 
          ProcessError::ProcessError(par_list), 
          Environ::Environ(par_list),
          aphid_name(aphid_name_) {};
    
    // using AphidPop::AphidPop;
    // using WaspPop::WaspPop;
    // using WaspAttack::WaspAttack;
    // using ProcessError::ProcessError;
    // using Environ::Environ;

    void show() const {
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

};









//' @export patch
class patch {
public:
    
    // -------
    // Members:
    // -------
    const uint harvest_period; // time points between harvests
    const uint harvest_offset; // time at which to begin harvests
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