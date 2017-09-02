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









// ====================================================================================
// ====================================================================================

// Classes for information that is constant through time

// ====================================================================================
// ====================================================================================




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
        arma::rowvec survs = arma::diagvec(leslie, -1).t();
        arma::rowvec fecs = leslie(0, arma::span(1, leslie.n_cols - 1));
        
        arma::uword N = std::min(static_cast<arma::uword>(6), leslie.n_rows);
        
        Rcout.precision(4);
        Rcout << std::fixed;
        
        Rcout << "< Aphid line constant info >" << endl;
        Rcout << "Resistances: (";
        Rcout << attack_surv(0) << ' ' << attack_surv(1) << ')' << endl;
        
        Rcout << "Survivals:   (";
        for (unsigned i = 0; i < N; i++) Rcout << survs(i) << ' ';
        if (leslie.n_rows > N) Rcout << "...";
        Rcout << ')' << endl;
        
        Rcout << "Fecundities: (";
        for (unsigned i = 0; i < N; i++) Rcout << fecs(i) << ' ';
        if (leslie.n_rows > N) Rcout << "...";
        Rcout << ')' << endl;
        
        
                
        return;
    }

};






#endif