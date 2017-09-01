# ifndef _ACERS_TYPES_H
# define _ACERS_TYPES_H

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

// // Aphid population information
// //' @export aphid_pop
// class aphid_pop {
// public:
//     // Fields
//     arma::mat leslie;       // Leslie matrix with survival and reproduction
//                             //   To create this, you need the following info:
//                             //     arma::uvec instar_days, double surv_juv, 
//                             //     arma::vec surv_adult, arma::vec repro
//     arma::vec X_0;          // initial aphid abundances by stage
//     double K;               // aphid density dependence
//     uint n_stages;          // number of aphid stages (i.e., days)
//     // Constructors
//     aphid_pop() {};
//     aphid_pop(double density_0, double K_, arma::uvec instar_days, 
//               double surv_juv, arma::vec surv_adult, arma::vec repro) {
//         leslie = leslie_matrix(instar_days, surv_juv, surv_adult, repro);
//         K = K_;
//         n_stages = arma::sum(instar_days);
//         X_0 = density_0 * leslie_sad(leslie);
//     };
//     
// };
// 
// 
// 
// 
// // Wasp population
// //' @export wasp_pop
// class wasp_pop {
// public:
//     // Fields
//     arma::vec Y_0;          // initial wasp abundances by stage
//     double sex_ratio;       // proportion of female wasps
//     double K_y;             // parasitized aphid density dependence
//     double s_y;             // parasitoid adult daily survival
//     arma::uvec mum_days;    // number of days per mummy stage (aphid alive & dead)
//     uint n_stages;          // number of wasp stages (i.e., days)
// 
//     // Constructors
//     wasp_pop() {};
//     wasp_pop(double density_0, double sex_ratio_, double K_y_,  double s_y_,
//              arma::uvec mum_days_) {
//         sex_ratio = sex_ratio_;
//         K_y = K_y_;
//         s_y = s_y_;
//         mum_days = mum_days_;
//         n_stages = arma::sum(mum_days_) + 1;
//         Y_0 = arma::zeros<arma::vec>(n_stages);
//         Y_0.tail(1) = density_0;
//     };
// 
// };
// 
// 
// 
// // Wasp attack
// //' @export wasp_attack
// class wasp_attack {
// public:
//     // Fields
//     arma::vec rel_attack;   // relative wasp attack rates by aphid stage
//     double a;               // overall parasitoid attack rate
//     double k;               // aggregation parameter of the nbinom distribution
//     double h;               // parasitoid attack rate handling time
//     arma::vec attack_surv;  // survival rates of singly & multiply attacked aphids
// 
//     // Constructors
//     wasp_attack() {};
//     wasp_attack(arma::vec rel_attack_, double a_, double k_, double h_,
//                 arma::vec attack_surv_) :
//         rel_attack(rel_attack_), a(a_), k(k_), h(h_), attack_surv(attack_surv_) {}
// };
// 
// 
// // Process error
// //' @export process_error
// class process_error {
// public:
//     // Fields
//     double sigma_x;         // environmental standard deviation for aphids
//     double sigma_y;         // environmental standard deviation for parasitoids
//     double rho;             // environmental correlation among instars
//     double demog_mult;      // Multiplier for demographic stochasticity
// 
//     // Constructors
//     process_error() {};
//     process_error(double sigma_x_, double sigma_y_, double rho_, double demog_mult_) :
//         sigma_x(sigma_x_), sigma_y(sigma_y_), rho(rho_), demog_mult(demog_mult_) {}
// };
// 
// 
// // Environment: Harvest, dispersal, and predation
// //' @export environ
// class environ {
// public:
//     // Fields
//     double harvest_surv;    // survival rate for living aphids during a harvest
//     double disp_aphid;      // dispersal rate for aphids
//     double disp_wasp;       // dispersal rate for wasps
//     arma::uvec disp_stages; // stages in which dispersal occurs in aphids
//     double pred_rate;       // predation on aphids and mummies
// 
//     // Constructors
//     environ() {};
//     environ(double harvest_surv_, double disp_aphid_, double disp_wasp_,
//             arma::uvec disp_stages_, double pred_rate_) :
//         harvest_surv(harvest_surv_), disp_aphid(disp_aphid_), disp_wasp(disp_wasp_),
//         disp_stages(disp_stages_), pred_rate(pred_rate_) {}
// };






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
    
    aphid_wasp_info(List par_list)
    {
        // Check that the list has all necessary parameters:
        CharacterVector needed_names = {
            "instar_days", "surv_juv", "surv_adult", "repro", "aphid_density_0",
            "K", "sex_ratio", "K_y", "s_y", "mum_days", "wasp_density_0",
            "rel_attack", "a", "k", "h", "attack_surv", "sigma_x", "sigma_y",
            "rho", "demog_mult", "harvest_surv", "disp_aphid", "disp_wasp",
            "disp_stages", "pred_rate"};
        CharacterVector list_names = par_list.names();
        for (auto n : needed_names) {
            if (std::find(list_names.begin(), list_names.end(), n) == list_names.end()) {
                Rcout << n << endl;
                stop("The above needed argument was not found");
            };
        }

        leslie = leslie_matrix(as<arma::uvec>(par_list["instar_days"]),
                               as<double>(par_list["surv_juv"]),
                               as<arma::vec>(par_list["surv_adult"]),
                               as<arma::vec>(par_list["repro"]));
        X_0 = as<double>(par_list["aphid_density_0"]) * leslie_sad(leslie);
        K = as<double>(par_list["K"]);
        n_aphid_stages = arma::sum(as<arma::uvec>(par_list["instar_days"]));

        sex_ratio = as<double>(par_list["sex_ratio"]);
        K_y = as<double>(par_list["K_y"]);
        s_y = as<double>(par_list["s_y"]);
        mum_days = as<arma::uvec>(par_list["mum_days"]);
        n_wasp_stages = arma::sum(mum_days) + 1;
        Y_0 = arma::zeros<arma::vec>(n_wasp_stages);
        Y_0.tail(1) = as<double>(par_list["wasp_density_0"]);;

        rel_attack = as<arma::vec>(par_list["rel_attack"]);
        a = as<double>(par_list["a"]);
        k = as<double>(par_list["k"]);
        h = as<double>(par_list["h"]);
        attack_surv = as<arma::vec>(par_list["attack_surv"]);

        sigma_x = as<double>(par_list["sigma_x"]);
        sigma_y = as<double>(par_list["sigma_y"]);
        rho = as<double>(par_list["rho"]);
        demog_mult = as<double>(par_list["demog_mult"]);

        harvest_surv = as<double>(par_list["harvest_surv"]);
        disp_aphid = as<double>(par_list["disp_aphid"]);
        disp_wasp = as<double>(par_list["disp_wasp"]);
        disp_stages = as<arma::uvec>(par_list["disp_stages"]);
        pred_rate = as<double>(par_list["pred_rate"]);
    }

    void show() {
        Rcout << "<< Aphid line constant info >>" << endl;
        Rcout << "Parasitoid resistance vector: ";
        attack_surv.t().print(Rcout);
        Rcout << "First 7x7 cells of the Leslie matrix:" << endl;
        leslie(arma::span(0,6),arma::span(0, 6)).print(Rcout);
        return;
    }

};






#endif