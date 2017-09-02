
Rcpp::loadModule("pop_nums_mod", TRUE)
Rcpp::loadModule("const_pop_mod", TRUE)
Rcpp::loadModule("aphid_line_mod", TRUE)


#' Create an const_pop object from input parameters.
#'
#' @param aphid_name Name of aphid line. Used for combining aphid numbers for dispersal
#'     among patches. Creating multiple `const_pop` objects with the `aphid_name`
#'     field but differing in other parameters simulates differences in environmental 
#'     effects.
#' @param aphid_density_0 Starting aphid density. Defaults to `populations$aphids_0`.
#' @param aphid_surv_juv Aphid juvenile survival. Defaults to `populations$surv_juv$high`.
#' @param aphid_surv_adult Vector of aphid adult survivals by stage. 
#'     Defaults to `populations$surv_adult$high`.
#' @param aphid_repro Vector of aphid reproductive rates by stage. 
#'     Defaults to `populations$repro$high`.
#' @param K Aphid density dependence. Defaults to `populations$K`.
#' @param wasp_density_0 Starting wasp density. Defaults to `populations$wasps_0`.
#' @param K_y Wasp density dependence. Defaults to `populations$K_y`.
#' @param s_y Wasp adult daily survival. Defaults to `populations$s_y`.
#' @param wasp_sex_ratio Proportion of female wasps. Defaults to `populations$sex_ratio`.
#' @param aphid_instar_days Vector of the number of stages (days) per aphid instar. 
#'     Defaults to `dev_times$instar_days$highT`.
#' @param mum_days Vector for the number of days per mummy stage. It should
#'     be of length 2, representing stages where the parasitized aphid is alive and 
#'     dead, respectively. Defaults to `dev_times$mum_days`.
#' @param rel_attack Vectof of relative wasp attack rates by aphid stage. 
#'     Defaults to `wasp_attack$rel_attack$highT`.
#' @param a Overall parasitoid attack rate. Defaults to `wasp_attack$a`.
#' @param k Aggregation parameter of the nbinom distribution. 
#'     Defaults to `wasp_attack$k`.
#' @param h Parasitoid attack rate handling time. 
#'     Defaults to `wasp_attack$h`.
#' @param attack_surv Vector of length 2 representing survival rates of singly & 
#'     multiply attacked aphids, respectively. Defaults to `numeric(2)`.
#' @param sigma_x Environmental standard deviation for aphids.
#'     Defaults to `environ$sigma_x`.
#' @param sigma_y Environmental standard deviation for parasitoids.
#'     Defaults to `environ$sigma_y`.
#' @param rho Environmental correlation among instars. Defaults to `environ$rho`.
#' @param demog_mult Multiplier for demographic stochasticity. Defaults to `1.0`.
#' @param harvest_surv Survival rate for living aphids during a harvest. 
#'     Defaults to `environ$harvest_surv`.
#' @param disp_aphid Dispersal rate for aphids. Defaults to `environ$disp_aphid`.
#' @param disp_wasp Dispersal rate for wasps. Defaults to `environ$disp_wasp`.
#' @param disp_start Stage in which dispersal starts in aphids. 
#'     Defaults to `environ$disp_start$high`.
#' @param pred_rate Predation on aphids and mummies. Defaults to `environ$pred_rate`.
#'
#' @return
#' 
#' Reference class 'Rcpp_const_pop' [package "sap"] with 24 fields
#' 
#' @export
#'
#' @examples
#' 
#' # Susceptible aphid line info
#' sus_line <- make_const_pop(
#'     "susceptible",
#'     aphid_density_0 = (1 - sap::populations$prop_resist) * 
#'         sap::populations$aphids_0)
#' sus_line
#' 
#' 
#' # Resistant aphid line info
#' # Differences from susceptible are resistance, lower reproduction, and 
#' # lower starting density
#' res_line <- make_const_pop(
#'     "resistant",
#'     attack_surv = sap::wasp_attack$attack_surv,
#'     aphid_density_0 = sap::populations$prop_resist * sap::populations$aphids_0,
#'     aphid_repro = sap::populations$repro$low)
#' res_line
#' 
make_const_pop <- function(
    aphid_name,
    aphid_density_0 = populations$aphids_0, 
    aphid_surv_juv = populations$surv_juv$high, 
    aphid_surv_adult = populations$surv_adult$high, 
    aphid_repro = populations$repro$high, 
    K = populations$K, 
    wasp_density_0 = populations$wasps_0, 
    K_y = populations$K_y, 
    s_y = populations$s_y, 
    wasp_sex_ratio = populations$sex_ratio, 
    aphid_instar_days = dev_times$instar_days$highT, 
    mum_days = dev_times$mum_days, 
    rel_attack = wasp_attack$rel_attack$highT, 
    a = wasp_attack$a, 
    k = wasp_attack$k, 
    h = wasp_attack$h, 
    attack_surv = numeric(2), 
    sigma_x = environ$sigma_x, 
    sigma_y = environ$sigma_y, 
    rho = environ$rho, 
    demog_mult = 1.0, 
    harvest_surv = environ$harvest_surv, 
    disp_aphid = environ$disp_aphid, 
    disp_wasp = environ$disp_wasp, 
    disp_start = environ$disp_start$highT,
    pred_rate = environ$pred_rate) {
    
    L = list()
    
    L[["aphid_name"]] = aphid_name;
    L[["instar_days"]] = aphid_instar_days;
    L[["surv_juv"]] = aphid_surv_juv;
    L[["surv_adult"]] = aphid_surv_adult;
    L[["repro"]] = aphid_repro;
    L[["aphid_density_0"]] = aphid_density_0;
    L[["K"]] = K;
    L[["wasp_density_0"]] = wasp_density_0;
    L[["K_y"]] = K_y;
    L[["s_y"]] = s_y;
    L[["sex_ratio"]] = wasp_sex_ratio;
    L[["mum_days"]] = mum_days;
    L[["rel_attack"]] = rel_attack;
    L[["a"]] = a;
    L[["k"]] = k;
    L[["h"]] = h;
    L[["attack_surv"]] = attack_surv;
    L[["sigma_x"]] = sigma_x;
    L[["sigma_y"]] = sigma_y;
    L[["rho"]] = rho;
    L[["demog_mult"]] = demog_mult;
    L[["harvest_surv"]] = harvest_surv;
    L[["disp_aphid"]] = disp_aphid;
    L[["disp_wasp"]] = disp_wasp;
    L[["disp_start"]] = disp_start;
    L[["pred_rate"]] = pred_rate;
    
    awi = new(const_pop, L)
    
    return(awi)
}


