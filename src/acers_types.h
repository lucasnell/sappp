# ifndef _ACERS_TYPES_H
# define _ACERS_TYPES_H

#include <RcppArmadillo.h> // arma namespace
#include <sitmo.h>         // parallel rng
#include <vector>          // vector class
#include <cmath>           // log, exp
#include <random>          // normal distribution
#include <cstdint>         // integer types
#include "math.h"          // leslie_matrix and leslie_sad


using namespace Rcpp;
using namespace std;
// using namespace arma;


typedef uint_fast8_t uint8;
typedef uint_fast32_t uint;
typedef int_fast32_t sint;
typedef uint_fast64_t uint64;
typedef int_fast64_t sint64;









// ====================================================================================
// ====================================================================================

// Classes for information that is constant through time

// ====================================================================================
// ====================================================================================

// Aphid population information
//' @export aphid_pop
class aphid_pop {
public:
    // Fields
    arma::mat leslie;       // Leslie matrix with survival and reproduction
                            //   To create this, you need the following info:
                            //     arma::uvec stage_days, double surv_juv, 
                            //     arma::vec surv_adult, arma::vec repro
    arma::vec X_0;          // initial aphid abundances by stage
    double K;               // aphid density dependence
    uint n_stages;          // number of aphid stages (i.e., days)
    // Constructors
    aphid_pop() {};
    aphid_pop(double density_0, double K_, arma::uvec stage_days, 
              double surv_juv, arma::vec surv_adult, arma::vec repro) {
        leslie = leslie_matrix(stage_days, surv_juv, surv_adult, repro);
        K = K_;
        n_stages = arma::sum(stage_days);
        X_0 = density_0 * leslie_sad(leslie);
    };
    
};




// Wasp population
//' @export wasp_pop
class wasp_pop {
public:
    // Fields
    arma::vec Y_0;          // initial wasp abundances by stage
    double sex_ratio;       // proportion of female wasps
    double K_y;             // parasitized aphid density dependence
    double s_y;             // parasitoid adult daily survival
    arma::uvec mum_days;    // number of days per mummy stage (aphid alive & dead)
    uint n_stages;          // number of wasp stages (i.e., days)

    // Constructors
    wasp_pop() {};
    wasp_pop(double density_0, double sex_ratio_, double K_y_,  double s_y_,
             arma::uvec mum_days_) {
        sex_ratio = sex_ratio_;
        K_y = K_y_;
        s_y = s_y_;
        mum_days = mum_days_;
        n_stages = arma::sum(mum_days_) + 1;
        Y_0 = arma::zeros<arma::vec>(n_stages);
        Y_0.tail(1) = density_0;
    };

};



// Wasp attack
//' @export wasp_attack
class wasp_attack {
public:
    // Fields
    arma::vec rel_attack;   // relative wasp attack rates by aphid stage
    double a;               // overall parasitoid attack rate
    double k;               // aggregation parameter of the nbinom distribution
    double h;               // parasitoid attack rate handling time
    arma::vec attack_surv;  // survival rates of singly & multiply attacked aphids

    // Constructors
    wasp_attack() {};
    wasp_attack(arma::vec rel_attack_, double a_, double k_, double h_,
                arma::vec attack_surv_) :
        rel_attack(rel_attack_), a(a_), k(k_), h(h_), attack_surv(attack_surv_) {}
};


// Process error
//' @export process_error
class process_error {
public:
    // Fields
    double sigma_x;         // environmental standard deviation for aphids
    double sigma_y;         // environmental standard deviation for parasitoids
    double rho;             // environmental correlation among instars
    double demog_mult;      // Multiplier for demographic stochasticity

    // Constructors
    process_error() {};
    process_error(double sigma_x_, double sigma_y_, double rho_, double demog_mult_) :
        sigma_x(sigma_x_), sigma_y(sigma_y_), rho(rho_), demog_mult(demog_mult_) {}
};


// Environment: Harvest, dispersal, and predation
//' @export environ
class environ {
public:
    // Fields
    double harvest_surv;    // survival rate for living aphids during a harvest
    double disp_aphid;      // dispersal rate for aphids
    double disp_wasp;       // dispersal rate for wasps
    arma::uvec disp_stages; // stages in which dispersal occurs in aphids
    double pred_rate;       // predation on aphids and mummies

    // Constructors
    environ() {};
    environ(double harvest_surv_, double disp_aphid_, double disp_wasp_,
            arma::uvec disp_stages_, double pred_rate_) :
        harvest_surv(harvest_surv_), disp_aphid(disp_aphid_), disp_wasp(disp_wasp_),
        disp_stages(disp_stages_), pred_rate(pred_rate_) {}
};






// Constant info about one aphid line and wasps that parasitize them
// This should be used among fields for the same aphid line

//' @export aphid_wasp_info
class aphid_wasp_info {
public:

    // --------
    // Parameters:
    // --------

    // Aphid population
    arma::mat leslie;       // Leslie matrix with survival and reproduction
    arma::vec X_0;          // initial aphid abundances by stage
    double K;               // aphid density dependence
    uint n_aphid_stages;    // number of aphid stages (i.e., days)

    // Wasp population
    arma::vec Y_0;          // initial wasp abundances by stage
    double sex_ratio;       // proportion of female wasps
    double K_y;             // parasitized aphid density dependence
    double s_y;             // parasitoid adult daily survival
    arma::uvec mum_days;    // number of days per mummy stage (aphid alive & dead)
    uint n_wasp_stages;     // number of wasp stages (i.e., days)

    // Wasp attack
    arma::vec rel_attack;   // relative wasp attack rates by aphid stage
    double a;               // overall parasitoid attack rate
    double k;               // aggregation parameter of the nbinom distribution
    double h;               // parasitoid attack rate handling time
    arma::vec attack_surv;  // survival rates of singly & multiply attacked aphids

    // Process error
    double sigma_x;         // environmental standard deviation for aphids
    double sigma_y;         // environmental standard deviation for parasitoids
    double rho;             // environmental correlation among instars
    double demog_mult;      // Multiplier for demographic stochasticity

    // Environment: Harvest, dispersal, and predation
    double harvest_surv;    // survival rate for living aphids during a harvest
    double disp_aphid;      // dispersal rate for aphids
    double disp_wasp;       // dispersal rate for wasps
    arma::uvec disp_stages; // stages in which dispersal occurs in aphids
    double pred_rate;       // predation on aphids and mummies


    // --------
    // Constructors:
    // --------

    aphid_wasp_info() {};
    
    aphid_wasp_info(aphid_pop aphid_pop_obj,
                    wasp_pop wasp_pop_obj,
                    wasp_attack wasp_attack_obj,
                    process_error process_error_obj,
                    environ environ_obj)
    {
        leslie = aphid_pop_obj.leslie;
        X_0 = aphid_pop_obj.X_0;
        K = aphid_pop_obj.K;
        n_aphid_stages = aphid_pop_obj.n_stages;

        Y_0 = wasp_pop_obj.Y_0;
        sex_ratio = wasp_pop_obj.sex_ratio;
        K_y = wasp_pop_obj.K_y;
        s_y = wasp_pop_obj.s_y;
        mum_days = wasp_pop_obj.mum_days;
        n_wasp_stages = wasp_pop_obj.n_stages;

        rel_attack = wasp_attack_obj.rel_attack;
        a = wasp_attack_obj.a;
        k = wasp_attack_obj.k;
        h = wasp_attack_obj.h;
        attack_surv = wasp_attack_obj.attack_surv;

        sigma_x = process_error_obj.sigma_x;
        sigma_y = process_error_obj.sigma_y;
        rho = process_error_obj.rho;
        demog_mult = process_error_obj.demog_mult;

        harvest_surv = environ_obj.harvest_surv;
        disp_aphid = environ_obj.disp_aphid;
        disp_wasp = environ_obj.disp_wasp;
        disp_stages = environ_obj.disp_stages;
        pred_rate = environ_obj.pred_rate;
    }

    void show() {
        Rcout << "<< Aphid line constant info >>" << endl;
        Rcout << "Parasitoid resistance vector:";
        for (auto s : attack_surv) Rcout << s;
        Rcout << endl;
        Rcout << "First 7x7 cells of the Leslie matrix:" << endl;
        leslie(arma::span(0,6),arma::span(0, 6)).print();
        return;
    }

};






#endif