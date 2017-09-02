
Rcpp::loadModule("pop_nums_mod", TRUE)
Rcpp::loadModule("const_pop_mod", TRUE)
Rcpp::loadModule("aphid_line_mod", TRUE)


# Create an const_pop object from input parameters
make_const_pop <- function(
    aphid_name,
    aphid_density_0, K, aphid_instar_days, aphid_surv_juv, aphid_surv_adult, 
    aphid_repro, wasp_density_0, K_y, s_y, mum_days, rel_attack, a, k, h, 
    sigma_x, sigma_y, rho, disp_stages,
    wasp_sex_ratio = 0.5, attack_surv = numeric(2), demog_mult = 1.0, 
    harvest_surv = 0.05, disp_aphid = 0.05, disp_wasp = 1.0, pred_rate = 0.8) {
    
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
    L[["disp_stages"]] = disp_stages;
    L[["pred_rate"]] = pred_rate;
    
    awi = new(const_pop, L)
    
    return(awi)
}


