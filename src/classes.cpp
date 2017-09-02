#include <RcppArmadillo.h> // arma namespace
#include <sitmo.h>         // parallel rng
#include <vector>          // vector class
#include <cmath>           // log, exp
#include <random>          // normal distribution
#include <cstdint>         // integer types

#include "sap_types.h"

using namespace Rcpp;
using namespace std;




RCPP_MODULE(sap_module) {
    
    class_<aphid_wasp_info>("aphid_wasp_info")
        .constructor<List>()
        .method("show", &aphid_wasp_info::show)
        .field_readonly("leslie", &aphid_wasp_info::leslie,
            "Leslie matrix with survival and reproduction")
        .field_readonly("X_0", &aphid_wasp_info::X_0, "initial aphid abundances by stage")
        .field_readonly("K", &aphid_wasp_info::K, "aphid density dependence")
        .field_readonly("n_aphid_stages", &aphid_wasp_info::n_aphid_stages,
            "number of aphid stages (i.e., days)")
        .field_readonly("Y_0", &aphid_wasp_info::Y_0, "initial wasp abundances by stage")
        .field_readonly("sex_ratio", &aphid_wasp_info::sex_ratio,
            "proportion of female wasps")
        .field_readonly("K_y", &aphid_wasp_info::K_y,
            "parasitized aphid density dependence")
        .field_readonly("s_y", &aphid_wasp_info::s_y, "parasitoid adult daily survival")
        .field_readonly("mum_days", &aphid_wasp_info::mum_days,
            "number of days per mummy stage (aphid alive & dead)")
        .field_readonly("n_wasp_stages", &aphid_wasp_info::n_wasp_stages,
            "number of wasp stages (i.e., days)")
        .field_readonly("rel_attack", &aphid_wasp_info::rel_attack,
            "relative wasp attack rates by aphid stage")
        .field_readonly("a", &aphid_wasp_info::a, "overall parasitoid attack rate")
        .field_readonly("k", &aphid_wasp_info::k,
            "aggregation parameter of the nbinom distribution")
        .field_readonly("h", &aphid_wasp_info::h, "parasitoid attack rate handling time")
        .field_readonly("attack_surv", &aphid_wasp_info::attack_surv,
            "survival rates of singly & multiply attacked aphids")
        .field_readonly("sigma_x", &aphid_wasp_info::sigma_x,
            "environmental standard deviation for aphids")
        .field_readonly("sigma_y", &aphid_wasp_info::sigma_y,
            "environmental standard deviation for parasitoids")
        .field_readonly("rho", &aphid_wasp_info::rho,
            "environmental correlation among instars")
        .field_readonly("demog_mult", &aphid_wasp_info::demog_mult,
            "Multiplier for demographic stochasticity")
        .field_readonly("harvest_surv", &aphid_wasp_info::harvest_surv,
            "survival rate for living aphids during a harvest")
        .field_readonly("disp_aphid", &aphid_wasp_info::disp_aphid,
            "dispersal rate for aphids")
        .field_readonly("disp_wasp", &aphid_wasp_info::disp_wasp,
            "dispersal rate for wasps")
        .field_readonly("disp_stages", &aphid_wasp_info::disp_stages,
            "stages in which dispersal occurs in aphids")
        .field_readonly("pred_rate", &aphid_wasp_info::pred_rate,
            "predation on aphids and mummies")
    ;
    
    
}
