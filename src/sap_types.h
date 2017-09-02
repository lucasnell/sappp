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










// Simple class to hold population numbers for aphids and wasps through time
//' @export pop_nums
class pop_nums {
public:

    // -------
    // Members
    // -------
    arma::vec X_t;              // Aphid density at time t
    arma::vec X_t1;             // Aphid density at time t+1
    arma::vec Y_t;              // Wasp density at time t
    arma::vec Y_t1;             // Wasp density at time t+1
    arma::vec A;                // Attack probabilities at time t
    
    // -------
    // Constructors
    // -------
    pop_nums(const arma::vec& X_0, const arma::vec& Y_0) 
        : X_t(X_0), X_t1(X_0), Y_t(Y_0), Y_t1(Y_0),
          A(arma::zeros<arma::vec>(X_0.n_elem)) {};
    
    pop_nums(const arma::uword& X_len, const arma::uword& Y_len) 
        : X_t(arma::zeros<arma::vec>(X_len)), 
          X_t1(arma::zeros<arma::vec>(X_len)), 
          Y_t(arma::zeros<arma::vec>(Y_len)), 
          Y_t1(arma::zeros<arma::vec>(Y_len)),
          A(arma::zeros<arma::vec>(X_len)) {};
    
    // -------
    // Methods
    // -------
    
    void show() const {
        
        Rcout << "Population numbers:" << endl;
        arma::uword N = std::min(static_cast<arma::uword>(6), X_t.n_elem);
        Rcout << "  * aphids[t]:     (";
        for (unsigned i = 0; i < N; i++) Rcout << X_t(i) << ' ';
        if (X_t.n_elem > N) Rcout << "...";
        Rcout << ')' << endl;
        Rcout << "  * aphids[t+1]:   (";
        for (unsigned i = 0; i < N; i++) Rcout << X_t1(i) << ' ';
        if (X_t1.n_elem > N) Rcout << "...";
        Rcout << ')' << endl;
        
        N = std::min(static_cast<arma::uword>(6), Y_t.n_elem);
        Rcout << "  * wasps[t]:      (";
        for (unsigned i = 0; i < N; i++) Rcout << Y_t(i) << ' ';
        if (Y_t.n_elem > N) Rcout << "...";
        Rcout << ')' << endl;
        Rcout << "  * wasps[t+1]:    (";
        for (unsigned i = 0; i < N; i++) Rcout << Y_t1(i) << ' ';
        if (Y_t1.n_elem > N) Rcout << "...";
        Rcout << ')' << endl;
        
        N = std::min(static_cast<arma::uword>(6), A.n_elem);
        Rcout << "  * Pr(attack)[t]: (";
        for (unsigned i = 0; i < N; i++) Rcout << A(i) << ' ';
        if (A.n_elem > N) Rcout << "...";
        Rcout << ')' << endl;
        
        return;
    }
    
};




// Constant info about one aphid line and wasps that parasitize them
// This should be used among fields for the same aphid line

//' @export const_pop
class const_pop {
public:

    // --------
    // Members:
    // --------

    string aphid_name;      // unique identifying name for this aphid line
    
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

    const_pop() {};
    
    const_pop(const List& par_list)
    {
        // Check that the list has all necessary parameters:
        CharacterVector needed_names = {
            "aphid_name",
            "instar_days", "surv_juv", "surv_adult", "repro", "aphid_density_0",
            "K", "sex_ratio", "K_y", "s_y", "mum_days", "wasp_density_0",
            "rel_attack", "a", "k", "h", "attack_surv", "sigma_x", "sigma_y",
            "rho", "demog_mult", "harvest_surv", "disp_aphid", "disp_wasp",
            "disp_stages", "pred_rate"};
        CharacterVector list_names = par_list.names();
        for (auto n : needed_names) {
            if (std::find(list_names.begin(), list_names.end(), n) == list_names.end()) {
                string error_msg = "The needed argument " + static_cast<string>(n) + 
                    " was not found";
                stop(error_msg);
            };
        }
        
        aphid_name = as<string>(par_list["aphid_name"]);

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

    void show() const {
        arma::rowvec survs = arma::diagvec(leslie, -1).t();
        arma::rowvec fecs = leslie(0, arma::span(1, leslie.n_cols - 1));
        
        arma::uword N = std::min(static_cast<arma::uword>(6), leslie.n_rows);
        
        Rcout.precision(4);
        Rcout << std::fixed;
        
        Rcout << "Constant info for '" << aphid_name << "' line" << endl;
        Rcout << "  * Resistances: (";
        Rcout << attack_surv(0) << ' ' << attack_surv(1) << ')' << endl;
        
        Rcout << "  * Survivals:   (";
        for (unsigned i = 0; i < N; i++) Rcout << survs(i) << ' ';
        if (leslie.n_rows > N) Rcout << "...";
        Rcout << ')' << endl;
        
        Rcout << "  * Fecundities: (";
        for (unsigned i = 0; i < N; i++) Rcout << fecs(i) << ' ';
        if (leslie.n_rows > N) Rcout << "...";
        Rcout << ')' << endl;
        
        
                
        return;
    }

};





// Info about an aphid line (and wasps that parasitize them), some of which changes 
// through time.
// This should be used for a single line in a single field.
// The pop_info field can use shared memory for all of the same line across fields

//' @export aphid_line
class aphid_line {
public:
    
    // --------
    // Members:
    // --------
    
    const const_pop& pop_info;  // Aphid line info (doesn't change through time)
    arma::vec X_t;              // Aphid density at time t
    arma::vec X_t1;             // Aphid density at time t+1
    arma::vec Y_t;              // Wasp density at time t
    arma::vec Y_t1;             // Wasp density at time t+1
    arma::vec A;                // Attack probabilities at time t+1
    
    
    // --------
    // Constructors:
    // --------
    
    aphid_line(const const_pop& info) 
        : pop_info(info),
          X_t(info.X_0), X_t1(info.X_0),
          Y_t(info.Y_0), Y_t1(info.Y_0),
          A(arma::zeros<arma::vec>(info.X_0.n_elem)) {};

    
    void show() const {
        
        Rcout << "< Aphid line info >" << endl;
        
        Rcout.precision(4);
        Rcout << std::fixed;
        arma::uword N;
        
        // pop_info.show(); // Do not do this!
        
        Rcout << endl;
        
        
        Rcout << "  Changing info:" << endl;
        
        N = std::min(static_cast<arma::uword>(6), X_t.n_elem);
        Rcout << "  * aphids[t]:     (";
        for (unsigned i = 0; i < N; i++) Rcout << X_t(i) << ' ';
        if (X_t.n_elem > N) Rcout << "...";
        Rcout << ')' << endl;
        Rcout << "  * aphids[t+1]:   (";
        for (unsigned i = 0; i < N; i++) Rcout << X_t1(i) << ' ';
        if (X_t1.n_elem > N) Rcout << "...";
        Rcout << ')' << endl;

        N = std::min(static_cast<arma::uword>(6), Y_t.n_elem);
        Rcout << "  * wasps[t]:      (";
        for (unsigned i = 0; i < N; i++) Rcout << Y_t(i) << ' ';
        if (Y_t.n_elem > N) Rcout << "...";
        Rcout << ')' << endl;
        Rcout << "  * wasps[t+1]:    (";
        for (unsigned i = 0; i < N; i++) Rcout << Y_t1(i) << ' ';
        if (Y_t1.n_elem > N) Rcout << "...";
        Rcout << ')' << endl;
        
        N = std::min(static_cast<arma::uword>(6), A.n_elem);
        Rcout << "  * Pr(attack)[t]: (";
        for (unsigned i = 0; i < N; i++) Rcout << A(i) << ' ';
        if (A.n_elem > N) Rcout << "...";
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
    arma::uvec harvest_times;  // times when harvesting occurs
    arma::uword max_time;      // number of time points to simulate
    double z;                  // Sum of all living aphids at time t
    double x;                  // Sum of non-parasitized aphids at time t
    double Y_m;                // Total number of adult wasps
    
    // -------
    // Constructor:
    // -------
    
    patch(const arma::uvec& harvest_times_, const arma::uword& max_time_)
        : harvest_times(harvest_times_), max_time(max_time_),
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