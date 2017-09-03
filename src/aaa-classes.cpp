#include <RcppArmadillo.h> // arma namespace
#include <sitmo.h>         // parallel rng
#include <vector>          // vector class
#include <cmath>           // log, exp
#include <random>          // normal distribution
#include <cstdint>         // integer types

#include "sap_types.h"

using namespace Rcpp;
using namespace std;

// RCPP_EXPOSED_CLASS(pop_nums)
// RCPP_MODULE(pop_nums_mod) {
//     
//     class_<pop_nums>("pop_nums")
//         .constructor<arma::vec, arma::vec>()
//         .method("show", &pop_nums::show)
//         .field("X_t", &pop_nums::X_t, "Aphid density at time t")
//         .field("X_t1", &pop_nums::X_t1, "Aphid density at time t+1")
//         .field("Y_t", &pop_nums::Y_t, "Wasp density at time t")
//         .field("Y_t1", &pop_nums::Y_t1, "Wasp density at time t+1")
//         .field("A", &pop_nums::A, "Attack probabilities at time t+1")
//     ;
// }

RCPP_EXPOSED_CLASS(const_pop)
RCPP_MODULE(const_pop_mod) {
    
    class_<const_pop>("const_pop")
        .constructor<List>()
        .method("show", &const_pop::show)
        .field_readonly("aphid_name", &const_pop::aphid_name, 
            "unique identifying name for this aphid line")
        // .field_readonly("leslie", &const_pop::leslie,
        //     "Leslie matrix with survival and reproduction")
        // .field_readonly("X_0", &const_pop::X_0, "initial aphid abundances by stage")
        // .field_readonly("K", &const_pop::K, "aphid density dependence")
        // .field_readonly("n_aphid_stages", &const_pop::n_aphid_stages,
        //     "number of aphid stages (i.e., days)")
        // .field_readonly("Y_0", &const_pop::Y_0, "initial wasp abundances by stage")
        // .field_readonly("sex_ratio", &const_pop::sex_ratio,
        //     "proportion of female wasps")
        // .field_readonly("K_y", &const_pop::K_y,
        //     "parasitized aphid density dependence")
        // .field_readonly("s_y", &const_pop::s_y, "parasitoid adult daily survival")
        // .field_readonly("mum_days", &const_pop::mum_days,
        //     "number of days per mummy stage (aphid alive & dead)")
        // .field_readonly("n_wasp_stages", &const_pop::n_wasp_stages,
        //     "number of wasp stages (i.e., days)")
        // .field_readonly("rel_attack", &const_pop::rel_attack,
        //     "relative wasp attack rates by aphid stage")
        // .field_readonly("a", &const_pop::a, "overall parasitoid attack rate")
        // .field_readonly("k", &const_pop::k,
        //     "aggregation parameter of the nbinom distribution")
        // .field_readonly("h", &const_pop::h, "parasitoid attack rate handling time")
        // .field_readonly("attack_surv", &const_pop::attack_surv,
        //     "survival rates of singly & multiply attacked aphids")
        // .field_readonly("sigma_x", &const_pop::sigma_x,
        //     "environmental standard deviation for aphids")
        // .field_readonly("sigma_y", &const_pop::sigma_y,
        //     "environmental standard deviation for parasitoids")
        // .field_readonly("rho", &const_pop::rho,
        //     "environmental correlation among instars")
        // .field_readonly("demog_mult", &const_pop::demog_mult,
        //     "Multiplier for demographic stochasticity")
        // .field_readonly("harvest_surv", &const_pop::harvest_surv,
        //     "survival rate for living aphids during a harvest")
        // .field_readonly("disp_aphid", &const_pop::disp_aphid,
        //     "dispersal rate for aphids")
        // .field_readonly("disp_wasp", &const_pop::disp_wasp,
        //     "dispersal rate for wasps")
        // .field_readonly("disp_start", &const_pop::disp_start,
        //     "stage in which dispersal starts in aphids")
        // .field_readonly("pred_rate", &const_pop::pred_rate,
        //     "predation on aphids and mummies")
    ;
}

// RCPP_MODULE(aphid_line_mod) {
//     
//     class_<aphid_line>("aphid_line")
//         .constructor<const_pop>()
//         .method("show", &aphid_line::show)
//         .field("X_t", &aphid_line::X_t, "Aphid density at time t")
//         .field("X_t1", &aphid_line::X_t1, "Aphid density at time t+1")
//         .field("Y_t", &aphid_line::Y_t, "Wasp density at time t")
//         .field("Y_t1", &aphid_line::Y_t1, "Wasp density at time t+1")
//         .field("A", &aphid_line::A, "Attack probabilities at time t+1")
//     ;
// }
