#include <RcppArmadillo.h> // arma namespace
#include <sitmo.h>         // parallel rng
#include <vector>          // vector class
#include <cmath>           // log, exp
#include <random>          // normal distribution
#include <cstdint>         // integer types

#include "acers_types.h"

using namespace Rcpp;
using namespace std;


// //[[Rcpp::export]]
// aphid_wasp_info make_aphid_wasp_info(
//         double aphid_density_0,
//         double K_,
//         arma::uvec aphid_stage_days,
//         double aphid_surv_juv,
//         arma::vec aphid_surv_adult,
//         arma::vec aphid_repro,
//         double wasp_density_0,
//         double K_y,
//         double s_y,
//         arma::uvec mum_days,
//         arma::vec rel_attack,
//         double a,
//         double k,
//         double h,
//         double sigma_x,
//         double sigma_y,
//         double rho,
//         arma::uvec disp_stages,
//         double wasp_sex_ratio = 0.5,
//         arma::vec attack_surv = NumericVector(2,0.0),
//         double demog_mult = 1.0,
//         double harvest_surv = 0.05,
//         double disp_aphid = 0.05,
//         double disp_wasp = 1.0,
//         double pred_rate = 0.8
//     ) {
// 
//     aphid_pop aphid_pop_obj(aphid_density_0, K_, aphid_stage_days, aphid_surv_juv,
//                             aphid_surv_adult, aphid_repro);
// 
//     wasp_pop wasp_pop_obj(wasp_density_0, wasp_sex_ratio, K_y,  s_y, mum_days);
// 
//     wasp_attack wasp_attack_obj(rel_attack, a, k, h, attack_surv);
// 
//     process_error process_error_obj(sigma_x, sigma_y, rho, demog_mult);
// 
//     environ environ_obj(harvest_surv, disp_aphid, disp_wasp, disp_stages, pred_rate);
// 
//     aphid_wasp_info aphid_wasp_info_obj(aphid_pop_obj, wasp_pop_obj, wasp_attack_obj,
//                                         process_error_obj, environ_obj);
// 
//     return aphid_wasp_info_obj;
// }








RCPP_MODULE(acers_module) {
    using namespace Rcpp;
    
    class_<aphid_pop>("aphid_pop")
        .constructor<double, double, arma::uvec, double, arma::vec, arma::vec>()
        .field("leslie", &aphid_pop::leslie, 
            "Leslie matrix with survival and reproduction")
        .field("X_0", &aphid_pop::X_0, "initial aphid abundances by stage")
        .field("K", &aphid_pop::K, "aphid density dependence")
        .field("n_stages", &aphid_pop::n_stages, "number of aphid stages (i.e., days)")
    ;
    class_<wasp_pop>( "wasp_pop")
        .constructor<double, double, double, double, arma::uvec>()
        .field("Y_0", &wasp_pop::Y_0, "initial wasp abundances by stage")
        .field("sex_ratio", &wasp_pop::sex_ratio, "proportion of female wasps")
        .field("K_y", &wasp_pop::K_y, "parasitized aphid density dependence")
        .field("s_y", &wasp_pop::s_y, "parasitoid adult daily survival")
        .field("mum_days", &wasp_pop::mum_days, 
            "number of days per mummy stage (aphid alive & dead)")
        .field("n_stages", &wasp_pop::n_stages, "number of wasp stages (i.e., days)")
    ;
    class_<wasp_attack>("wasp_attack")
        .constructor<arma::vec, double, double, double, arma::vec>()
        .field("rel_attack", &wasp_attack::rel_attack, 
            "relative wasp attack rates by aphid stage")
        .field("a", &wasp_attack::a, "overall parasitoid attack rate")
        .field("k", &wasp_attack::k, "aggregation parameter of the nbinom distribution")
        .field("h", &wasp_attack::h, "parasitoid attack rate handling time")
        .field("attack_surv", &wasp_attack::attack_surv, 
            "survival rates of singly & multiply attacked aphids")
    
    ;
    class_<process_error>("process_error")
        .constructor<double, double, double, double>()
    
        .field("sigma_x", &process_error::sigma_x, 
            "environmental standard deviation for aphids")
        .field("sigma_y", &process_error::sigma_y, 
            "environmental standard deviation for parasitoids")
        .field("rho", &process_error::rho, "environmental correlation among instars")
        .field("demog_mult", &process_error::demog_mult, 
            "Multiplier for demographic stochasticity")
    ;
    class_<environ>( "environ")
        .constructor<double, double, double, arma::uvec, double>()
        .field("harvest_surv", &environ::harvest_surv, 
            "survival rate for living aphids during a harvest")
        .field("disp_aphid", &environ::disp_aphid, "dispersal rate for aphids")
        .field("disp_wasp", &environ::disp_wasp, "dispersal rate for wasps")
        .field("disp_stages", &environ::disp_stages, 
            "stages in which dispersal occurs in aphids")
        .field("pred_rate", &environ::pred_rate, "predation on aphids and mummies")
    ;
    class_<aphid_wasp_info>("aphid_wasp_info")
        .constructor<aphid_pop, wasp_pop, wasp_attack, process_error, environ>()
    //     .method("show", &aphid_wasp_info::show)
    //     .field_readonly("leslie", &aphid_wasp_info::leslie,
    //         "Leslie matrix with survival and reproduction")
    //     .field_readonly("X_0", &aphid_wasp_info::X_0, "initial aphid abundances by stage")
    //     .field_readonly("K", &aphid_wasp_info::K, "aphid density dependence")
    //     .field_readonly("n_aphid_stages", &aphid_wasp_info::n_aphid_stages,
    //         "number of aphid stages (i.e., days)")
    //     .field_readonly("Y_0", &aphid_wasp_info::Y_0, "initial wasp abundances by stage")
    //     .field_readonly("sex_ratio", &aphid_wasp_info::sex_ratio,
    //         "proportion of female wasps")
    //     .field_readonly("K_y", &aphid_wasp_info::K_y,
    //         "parasitized aphid density dependence")
    //     .field_readonly("s_y", &aphid_wasp_info::s_y, "parasitoid adult daily survival")
    //     .field_readonly("mum_days", &aphid_wasp_info::mum_days,
    //         "number of days per mummy stage (aphid alive & dead)")
    //     .field_readonly("n_wasp_stages", &aphid_wasp_info::n_wasp_stages,
    //         "number of wasp stages (i.e., days)")
    //     .field_readonly("rel_attack", &aphid_wasp_info::rel_attack,
    //         "relative wasp attack rates by aphid stage")
    //     .field_readonly("a", &aphid_wasp_info::a, "overall parasitoid attack rate")
    //     .field_readonly("k", &aphid_wasp_info::k,
    //         "aggregation parameter of the nbinom distribution")
    //     .field_readonly("h", &aphid_wasp_info::h, "parasitoid attack rate handling time")
    //     .field_readonly("attack_surv", &aphid_wasp_info::attack_surv,
    //         "survival rates of singly & multiply attacked aphids")
    //     .field_readonly("sigma_x", &aphid_wasp_info::sigma_x,
    //         "environmental standard deviation for aphids")
    //     .field_readonly("sigma_y", &aphid_wasp_info::sigma_y,
    //         "environmental standard deviation for parasitoids")
    //     .field_readonly("rho", &aphid_wasp_info::rho,
    //         "environmental correlation among instars")
    //     .field_readonly("demog_mult", &aphid_wasp_info::demog_mult,
    //         "Multiplier for demographic stochasticity")
    //     .field_readonly("harvest_surv", &aphid_wasp_info::harvest_surv,
    //         "survival rate for living aphids during a harvest")
    //     .field_readonly("disp_aphid", &aphid_wasp_info::disp_aphid,
    //         "dispersal rate for aphids")
    //     .field_readonly("disp_wasp", &aphid_wasp_info::disp_wasp,
    //         "dispersal rate for wasps")
    //     .field_readonly("disp_stages", &aphid_wasp_info::disp_stages,
    //         "stages in which dispersal occurs in aphids")
    //     .field_readonly("pred_rate", &aphid_wasp_info::pred_rate,
    //         "predation on aphids and mummies")
    ;
    
    
}
