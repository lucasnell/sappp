
Rcpp::loadModule("AphidWasp_module", TRUE)


#' Create an AphidWasp object from input parameters.
#' 
#' Most arguments are \code{NULL} by defaults, which means the output object
#' will have values defined in \code{vignette("parameters")}.
#'
#' @param aphid_name Name of aphid line. Used for combining aphid numbers for dispersal
#'     among patches. Creating multiple `AphidWasp` objects with the `aphid_name`
#'     field but differing in other parameters simulates differences in environmental 
#'     effects.
#' @param aphid_density_0 Starting aphid density. 
#'     Defaults to \code{NULL}.
#' @param aphid_surv_juv Aphid juvenile survival. 
#'     Defaults to \code{NULL}.
#' @param aphid_surv_adult Vector of aphid adult survivals by stage. 
#'     Defaults to \code{NULL}.
#' @param aphid_repro Vector of aphid reproductive rates by stage. 
#'     Defaults to \code{NULL}.
#' @param K Aphid density dependence. Defaults to \code{NULL}.
#' @param wasp_density_0 Starting wasp density. Defaults to \code{NULL}.
#' @param K_y Wasp density dependence. Defaults to \code{NULL}.
#' @param s_y Wasp adult daily survival. Defaults to \code{NULL}.
#' @param wasp_sex_ratio Proportion of female wasps. 
#'     Defaults to \code{NULL}.
#' @param aphid_instar_days Vector of the number of stages (days) per aphid instar. 
#'     Defaults to \code{NULL}.
#' @param mum_days Vector for the number of days per mummy stage. It should
#'     be of length 2, representing stages where the parasitized aphid is alive and 
#'     dead, respectively. Defaults to \code{NULL}.
#' @param rel_attack Vectof of relative wasp attack rates by aphid stage. 
#'     Defaults to \code{NULL}.
#' @param a Overall parasitoid attack rate. Defaults to \code{NULL}.
#' @param k Aggregation parameter of the nbinom distribution. 
#'     Defaults to \code{NULL}.
#' @param h Parasitoid attack rate handling time. 
#'     Defaults to \code{NULL}.
#' @param attack_surv Vector of length 2 representing survival rates of singly & 
#'     multiply attacked aphids, respectively. Defaults to `numeric(2)`.
#' @param sigma_x Environmental standard deviation for aphids.
#'     Defaults to \code{NULL}.
#' @param sigma_y Environmental standard deviation for parasitoids.
#'     Defaults to \code{NULL}.
#' @param rho Environmental correlation among instars. Defaults to \code{NULL}.
#' @param demog_mult Multiplier for demographic stochasticity. Defaults to `1.0`.
#' @param harvest_surv Survival rate for living aphids during a harvest. 
#'     Defaults to \code{NULL}.
#' @param disp_aphid Dispersal rate for aphids. Defaults to \code{NULL}.
#' @param disp_wasp Dispersal rate for wasps. Defaults to \code{NULL}.
#' @param disp_start Stage in which dispersal starts in aphids. 
#'     Defaults to \code{NULL}.
#' @param pred_rate Predation on aphids and mummies. Defaults to \code{NULL}.
#'
#' @return
#' 
#' \code{Reference class 'Rcpp_AphidWasp' [package "sap"] with 29 fields}
#' 
#' @export
#'
#' @examples
#' 
#' # Susceptible aphid line info
#' sus_line <- make_AphidWasp(
#'     "susceptible",
#'     aphid_density_0 = (1 - sap::populations$prop_resist) * 
#'         sap::populations$aphids_0)
#' sus_line
#' 
#' 
#' # Resistant aphid line info
#' # Differences from susceptible are resistance, lower reproduction, and 
#' # lower starting density
#' res_line <- make_AphidWasp(
#'     "resistant",
#'     attack_surv = sap::wasp_attack$attack_surv,
#'     aphid_density_0 = sap::populations$prop_resist * sap::populations$aphids_0,
#'     aphid_repro = sap::populations$repro$low)
#' res_line
#' 
make_AphidWasp <- function(
    aphid_name,
    aphid_density_0 = NULL, 
    aphid_surv_juv = NULL, 
    aphid_surv_adult = NULL, 
    aphid_repro = NULL, 
    K = NULL, 
    wasp_density_0 = NULL, 
    K_y = NULL, 
    s_y = NULL, 
    wasp_sex_ratio = NULL, 
    aphid_instar_days = NULL, 
    mum_days = NULL, 
    rel_attack = NULL, 
    a = NULL, 
    k = NULL, 
    h = NULL, 
    attack_surv = numeric(2), 
    sigma_x = NULL, 
    sigma_y = NULL, 
    rho = NULL, 
    demog_mult = 1.0, 
    harvest_surv = NULL, 
    disp_aphid = NULL, 
    disp_wasp = NULL, 
    disp_start = NULL, 
    pred_rate = NULL) {
    
    if (is.null(aphid_density_0)) aphid_density_0 = sap::populations$aphids_0
    if (is.null(aphid_surv_juv)) aphid_surv_juv = sap::populations$surv_juv$high
    if (is.null(aphid_surv_adult)) aphid_surv_adult = sap::populations$surv_adult$high
    if (is.null(aphid_repro)) aphid_repro = sap::populations$repro$high
    if (is.null(K)) K = sap::populations$K
    if (is.null(wasp_density_0)) wasp_density_0 = sap::populations$wasps_0
    if (is.null(K_y)) K_y = sap::populations$K_y
    if (is.null(s_y)) s_y = sap::populations$s_y
    if (is.null(wasp_sex_ratio)) wasp_sex_ratio = sap::populations$sex_ratio
    if (is.null(aphid_instar_days)) aphid_instar_days = sap::dev_times$instar_days$highT
    if (is.null(mum_days)) mum_days = sap::dev_times$mum_days
    if (is.null(rel_attack)) rel_attack = sap::wasp_attack$rel_attack$highT
    if (is.null(a)) a = sap::wasp_attack$a
    if (is.null(k)) k = sap::wasp_attack$k
    if (is.null(h)) h = sap::wasp_attack$h
    if (is.null(sigma_x)) sigma_x = sap::environ$sigma_x
    if (is.null(sigma_y)) sigma_y = sap::environ$sigma_y
    if (is.null(rho)) rho = sap::environ$rho
    if (is.null(harvest_surv)) harvest_surv = sap::environ$harvest_surv
    if (is.null(disp_aphid)) disp_aphid = sap::environ$disp_aphid
    if (is.null(disp_wasp)) disp_wasp = sap::environ$disp_wasp
    if (is.null(disp_start)) disp_start = sap::environ$disp_start$highT
    if (is.null(pred_rate)) pred_rate = sap::environ$pred_rate
    
    L = list()
    
    # L[["aphid_name"]] = aphid_name;
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
    
    awi = new(AphidWasp, aphid_name, L)
    
    return(awi)
}


