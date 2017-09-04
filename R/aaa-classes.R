
Rcpp::loadModule("sap_module", TRUE)


#' Create a list for creating an AphidWasp object.
#' 
#' Most arguments are \code{NULL} by default, which means the output object
#' will have values defined in \code{vignette("parameters")}.
#' When multiple values are present in the \code{vignette("parameters")} parameters,
#' defaults here are those for high growth rates and high temperatures.
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
#'     multiply attacked aphids, respectively. Defaults to \code{numeric(2)}.
#' @param sigma_x Environmental standard deviation for aphids.
#'     Defaults to \code{NULL}.
#' @param sigma_y Environmental standard deviation for parasitoids.
#'     Defaults to \code{NULL}.
#' @param rho Environmental correlation among instars. Defaults to \code{NULL}.
#' @param demog_mult Multiplier for demographic stochasticity. Defaults to \code{1.0}.
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
#' List with all necessary info to create an AphidWasp object
#' 
#'
#' 
AphidWasp_list <- function(
    aphid_name, aphid_density_0 = NULL, aphid_surv_juv = NULL, aphid_surv_adult = NULL, 
    aphid_repro = NULL, K = NULL, wasp_density_0 = NULL, K_y = NULL, s_y = NULL, 
    wasp_sex_ratio = NULL, aphid_instar_days = NULL, mum_days = NULL, rel_attack = NULL, 
    a = NULL, k = NULL, h = NULL, attack_surv = numeric(2), sigma_x = NULL, 
    sigma_y = NULL, rho = NULL, demog_mult = 1.0, harvest_surv = NULL, disp_aphid = NULL,
    disp_wasp = NULL, disp_start = NULL, pred_rate = NULL) {
    
    if (is.null(aphid_density_0)) aphid_density_0 <- sap::populations$aphids_0
    if (is.null(aphid_surv_juv)) aphid_surv_juv <- sap::populations$surv_juv$high
    if (is.null(aphid_surv_adult)) aphid_surv_adult <- sap::populations$surv_adult$high
    if (is.null(aphid_repro)) aphid_repro <- sap::populations$repro$high
    if (is.null(K)) K <- sap::populations$K
    if (is.null(wasp_density_0)) wasp_density_0 <- sap::populations$wasps_0
    if (is.null(K_y)) K_y <- sap::populations$K_y
    if (is.null(s_y)) s_y <- sap::populations$s_y
    if (is.null(wasp_sex_ratio)) wasp_sex_ratio <- sap::populations$sex_ratio
    if (is.null(aphid_instar_days)) aphid_instar_days <- sap::dev_times$instar_days$highT
    if (is.null(mum_days)) mum_days <- sap::dev_times$mum_days
    if (is.null(rel_attack)) rel_attack <- sap::wasp_attack$rel_attack$highT
    if (is.null(a)) a <- sap::wasp_attack$a
    if (is.null(k)) k <- sap::wasp_attack$k
    if (is.null(h)) h <- sap::wasp_attack$h
    if (is.null(sigma_x)) sigma_x <- sap::environ$sigma_x
    if (is.null(sigma_y)) sigma_y <- sap::environ$sigma_y
    if (is.null(rho)) rho <- sap::environ$rho
    if (is.null(harvest_surv)) harvest_surv <- sap::environ$harvest_surv
    if (is.null(disp_aphid)) disp_aphid <- sap::environ$disp_aphid
    if (is.null(disp_wasp)) disp_wasp <- sap::environ$disp_wasp
    if (is.null(disp_start)) disp_start <- sap::environ$disp_start$highT
    if (is.null(pred_rate)) pred_rate <- sap::environ$pred_rate
    
    
    if (is.character(aphid_surv_juv)) {
        stopifnot(aphid_surv_juv %in% names(sap::populations$surv_juv))
        aphid_surv_juv <- sap::populations$surv_juv[[aphid_surv_juv]]
    }
    if (is.character(aphid_surv_adult)) {
        stopifnot(aphid_surv_adult %in% names(sap::populations$surv_adult))
        aphid_surv_adult <- sap::populations$surv_adult[[aphid_surv_adult]]
    }
    if (is.character(aphid_repro)) {
        stopifnot(aphid_repro %in% names(sap::populations$repro))
        aphid_repro <- sap::populations$repro[[aphid_repro]]
    }
    if (is.character(attack_surv)) {
        stopifnot(attack_surv %in% c('resistant', 'susceptible'))
        if (attack_surv == 'resistant') {
            attack_surv <- wasp_attack$attack_surv
        } else if (attack_surv == 'susceptible') {
            attack_surv <- c(0.0, 0.0)
        }
    } else if (is.numeric(attack_surv) & length(attack_surv) == 1) {
        stopifnot(attack_surv >= 0)
        attack_surv <- wasp_attack$attack_surv * attack_surv
        attack_surv <- ifelse(attack_surv > 1, 1, attack_surv)
    }
    
    
    # Rcpp modules (C++ classes Rcpp converts to R reference classes)
    # can only take <= 6 input arguments, so I'm combining all these arguments into a 
    # single list
    L = list()
    L[["aphid_name"]] = aphid_name
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
    
    
    # awi = new(AphidWasp, L)
    
    return(L)
}







# Create all parameter lists for all patches and one aphid line
one_pop_lists <- function(.aphid_name, .n_patches, ...) {
    if (missing(...)) {
        return(rep(list(list(aphid_name = .aphid_name)), .n_patches))
    }
    # Repeat each argument .n_patches times
    par_list <- lapply(
        list(...), 
        function(.par) {
            if (is.list(.par)) {
                if ((.n_patches %% length(.par)) != 0 | length(.par) == 0) {
                    stop(paste("length of all arguments must be a non-zero",
                               "factor of .n_patches"))
                }
                n_reps <- .n_patches / length(.par)
                outl <- lapply(.par, function(y) rep(list(y), n_reps))
                outl <- unlist(outl, recursive = FALSE)
            } else {
                outl <- rep(list(.par), .n_patches)
            }
            return(outl)
        })
    # Combine each argument so I now have a list with .n_patches items, each of 
    # which has all the input arguments
    par_list <- lapply(1:.n_patches, function(i) {
        outl <- c(.aphid_name, lapply(names(par_list), function(n) par_list[[n]][[i]]))
        names(outl) <- c("aphid_name", names(par_list))
        return(outl)
    })
    
    par_lists <- lapply(par_list, do.call, what = AphidWasp_list)
    
    return(par_lists)
}



# Create all parameter lists for all patches and one aphid line
# Each named entry must correspond with that from AphidWasp_list (e.g., `K = 0.5`).
# If you want different values per aphid-wasp combo, you should put the different 
# values in a list (e.g., `K = list(0.5, 0.6)`).
# If you want it to vary among and within aphid-wasp combos for >= 1 combo 
# (e.g., to simulate environmental effects), nest the list another level.
# For example, `K = list(0.5, list(0.55, 0.60))` would vary by combo, but would also
# vary among patches in the second aphid-wasp combo.
# The number of values among all aphid-wasp combos must be a factor of the number of 
# combos.
# Similarly, the number of values among all a combo's patches must be a factor of the
# number of patches. The number of patches is the same among combos.
all_pop_lists <- function(.n_pops, .n_patches, ..., .aphid_names = NULL) {
    
    stopifnot(.n_pops > 0, .n_patches > 0)
    
    if (is.null(.aphid_names)) {
        .aphid_names <- sprintf(paste0("pop_%0", ceiling(log10(.n_pops + 1)), "d"), 
                                1:.n_pops)
    } else {
        stopifnot(length(.aphid_names) == .n_pops)
    }
    
    if (missing(...)) {
        lapply(.aphid_names, function(.n) multiply_par_lists(.n, .n_patches))
    }
    # Repeat each argument .n_pops times
    par_list <- lapply(
        list(...),
        function(.par) {
            if (is.list(.par)) {
                if ((.n_pops %% length(.par)) != 0 | length(.par) == 0) {
                    stop(paste("length of all arguments must be a non-zero",
                               "factor of .n_pops"))
                }
                n_reps <- .n_pops / length(.par)
                outl <- lapply(.par, function(y) rep(list(y), n_reps))
                outl <- unlist(outl, recursive = FALSE)
            } else {
                outl <- rep(list(.par), .n_pops)
            }
            return(outl)
        })
    # Combine each argument so I now have a list with .n_pops items, each of 
    # which has all the input arguments
    par_list <- lapply(1:.n_pops, function(i) {
        outl <- c(.aphid_names[i], .n_patches, 
                  lapply(names(par_list), function(n) par_list[[n]][[i]]))
        names(outl) <- c(".aphid_name", ".n_patches", names(par_list))
        return(outl)
    })
    
    par_lists <- lapply(par_list, do.call, what = one_pop_lists)
    
    # Now reorganize it so that the list is nested by patch, then by pop
    par_lists <- lapply(1:.n_patches, 
                        function(pati) {
                            lapply(1:.n_pops, 
                                   function(popi) {
                                       par_lists[[popi]][[pati]]
                                   })
                        })
    
    return(par_lists)
}




