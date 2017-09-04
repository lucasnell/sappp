#include <RcppArmadillo.h> // arma namespace
#include <sitmo.h>         // parallel rng
#include <vector>          // vector class
#include <cmath>           // log, exp
#include <random>          // normal distribution
#include <cstdint>         // integer types

#include "sap_types.h"

using namespace Rcpp;
using namespace std;


RCPP_EXPOSED_CLASS(AphidWasp)
RCPP_EXPOSED_CLASS(SimSummary)
RCPP_EXPOSED_CLASS(OnePatch)
RCPP_EXPOSED_CLASS(SimPatches)

    



RCPP_MODULE(sap_module) {
    
    class_<SimSummary>("SimSummary")
        .field("aphids", &SimSummary::aphids, "density of unparasitized aphids")
        .field("parasit", &SimSummary::parasit, 
            "density of parasitized, but alive, aphids")
        .field("mummies", &SimSummary::mummies, "density of mummies")
        .field("wasps", &SimSummary::wasps, "density of adult wasps")
    ;
    
    class_<AphidPop>("AphidPop")
       .field_readonly("leslie", &AphidPop::leslie,
            "Leslie matrix with survival and reproduction")
        .field_readonly("X_0", &AphidPop::X_0, "initial aphid abundances by stage")
        .field_readonly("K", &AphidPop::K, "aphid density dependence")
        .field_readonly("n_stages", &AphidPop::n_aphid_stages,
            "number of aphid stages (i.e., days)")
        .field("X_t", &AphidPop::X_t, "Aphid density at time t")
        .field("X_t1", &AphidPop::X_t1, "Aphid density at time t+1")
    ;
    class_<WaspPop>("WaspPop")
        .field_readonly("Y_0", &WaspPop::Y_0, "initial wasp abundances by stage")
        .field_readonly("sex_ratio", &WaspPop::sex_ratio, "proportion of female wasps")
        .field_readonly("K_y", &WaspPop::K_y, "parasitized aphid density dependence")
        .field_readonly("s_y", &WaspPop::s_y, "parasitoid adult daily survival")
        .field_readonly("mum_days", &WaspPop::mum_days,
            "number of days per mummy stage (aphid alive & dead)")
        .field_readonly("n_stages", &WaspPop::n_wasp_stages,
            "number of wasp stages (i.e., days)")
        .field("Y_t", &WaspPop::Y_t, "Wasp density at time t")
        .field("Y_t1", &WaspPop::Y_t1, "Wasp density at time t+1")
    ;
    
    class_<WaspAttack>("WaspAttack")
        .field_readonly("rel_attack", &WaspAttack::rel_attack,
            "relative wasp attack rates by aphid stage")
        .field_readonly("a", &WaspAttack::a, "overall parasitoid attack rate")
        .field_readonly("k", &WaspAttack::k,
            "aggregation parameter of the nbinom distribution")
        .field_readonly("h", &WaspAttack::h, "parasitoid attack rate handling time")
        .field_readonly("attack_surv", &WaspAttack::attack_surv,
            "survival rates of singly & multiply attacked aphids")
        .field("A", &WaspAttack::A, "attack probabilities at time t")
    ;
    
    class_<ProcessError>("ProcessError")
        .field_readonly("sigma_x", &ProcessError::sigma_x,
            "environmental standard deviation for aphids")
        .field_readonly("sigma_y", &ProcessError::sigma_y,
            "environmental standard deviation for parasitoids")
        .field_readonly("rho", &ProcessError::rho,
            "environmental correlation among instars")
        .field_readonly("demog_mult", &ProcessError::demog_mult,
            "multiplier for demographic stochasticity")
        ;
    
    class_<Environ>("Environ")
        .field_readonly("harvest_surv", &Environ::harvest_surv,
            "survival rate for living aphids during a harvest")
        .field_readonly("disp_aphid", &Environ::disp_aphid, "dispersal rate for aphids")
        .field_readonly("disp_wasp", &Environ::disp_wasp, "dispersal rate for wasps")
        .field_readonly("disp_start", &Environ::disp_start,
            "stage in which dispersal starts in aphids")
        .field_readonly("pred_rate", &Environ::pred_rate,
            "predation on aphids and mummies")
        ;
    
    
    class_<AphidWasp>("AphidWasp")
        .derives<AphidPop>("AphidPop")
        .derives<WaspPop>("WaspPop")
        .derives<WaspAttack>("WaspAttack")
        .derives<ProcessError>("ProcessError")
        .derives<Environ>("Environ")
        .constructor<List>()
        .method("show", &AphidWasp::show)
        .field_readonly("aphid_name", &AphidWasp::aphid_name, 
            "unique identifying name for this aphid line")
    ;
    
    class_<OnePatch>("OnePatch")
        .field_readonly("harvest_period", &OnePatch::harvest_period, 
            "time points between harvests")
        .field_readonly("harvest_offset", &OnePatch::harvest_offset, 
            "time at which to begin harvests")
        .field("z", &OnePatch::z, "Sum of all living aphids at time t")
        .field("x", &OnePatch::x, "Sum of non-parasitized aphids at time t")
        .field("Y_m", &OnePatch::Y_m, "Total number of adult wasps")
        .method("show", &OnePatch::show)
    ;
    
    class_<SimPatches>("SimPatches")
        .constructor<vector<vector<List>>,vector<uint>,vector<uint>>()
        .field_readonly("aphid_names", &SimPatches::aphid_names, 
            "vector of aphid names (same for all patches)")
        .field("t", &SimPatches::t, "current time")
        .method("simulate", &SimPatches::simulate,
            "reset object and simulate a set number of time periods")
    ;
}

