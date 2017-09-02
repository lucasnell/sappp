
Rcpp::loadModule("pop_nums_mod", TRUE)
Rcpp::loadModule("const_pop_mod", TRUE)
Rcpp::loadModule("aphid_line_mod", TRUE)


# Create an const_pop object from input parameters
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
    disp_start = environ$disp_start$highT,
    harvest_surv = environ$harvest_surv, 
    disp_aphid = environ$disp_aphid, 
    disp_wasp = environ$disp_wasp, 
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


