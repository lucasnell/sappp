

# See \code{vignette("variability")} for how this was done
# For .mean that have values on the bounds (e.g., 0 or 1 for "aphid_surv_adult"), 
# I'm assuming those zeros are absolute. 
# For example, we can be quite sure that aphids do not live past 30 days, 
# so I don't want to add variability to this.

make_n_values <- function(.par_name, .n_pops, .n_patches, .sd_pops, 
                          .mean = NULL, .sd_patches = 0) {
    par_programmed <- c('aphid_density_0', 'aphid_surv_juv', 'aphid_surv_adult', 
                        'aphid_repro', 'K', 'disp_aphid', 'attack_surv')
    if (!.par_name %in% par_programmed) {
        stop(paste("parameter name", .par_name, "isn't programmed to be varied.",
                   "You have to do it manually."))
    }
    
    default_means <- list(
        aphid_density_0 = sappp::populations$aphids_0,
        aphid_surv_juv = Reduce(`+`, sappp::populations$surv_juv)/2,
        aphid_surv_adult = Reduce(`+`, sappp::populations$surv_adult)/2,
        aphid_repro = Reduce(`+`, sappp::populations$repro)/2,
        K = sappp::populations$K,
        disp_aphid = sappp::environ$disp_aphid,
        attack_surv = sappp::wasp_attack$attack_surv / 2
    )
    if (is.null(.mean)) .mean <- default_means[[.par_name]]
    # One copy is sufficient if you don't want any variability
    if (.sd_pops == 0 & .sd_patches == 0) return(.mean)
    
    # Parameters bound by (0,Inf)
    if (.par_name %in% c('aphid_density_0', 'aphid_repro')) {
        trans <- log
        inv_trans <- exp
    # Parameters bound by (0,1)
    } else {
        trans <- logit
        inv_trans <- inv_logit
    }
    
    # One random deviate value
    one_dev <- function(tm, .sd) {
        out <- inv_trans(tm + rnorm(1, sd = .sd))
        dim(out) <- dim(.mean)
        return(out)
    }
    
    .trans_mean <- trans(.mean)
    # Random deviates
    n_values <- replicate(.n_pops, one_dev(.trans_mean, .sd_pops), simplify = FALSE)
    
    if (.sd_patches > 0) {
        n_values <- lapply(n_values, 
                           function(x) {
                               replicate(.n_patches, one_dev(trans(x), .sd_patches), 
                                         simplify = FALSE)
                           })
    }
    
    return(n_values)
    
}







#' Create SimPatches object.
#' 
#' For parameters \code{aphid_density_0}, \code{aphid_surv_juv}, \code{aphid_surv_adult},
#' \code{aphid_repro}, \code{K}, \code{disp_aphid}, and \code{attack_surv}, 
#' you should provide a list containing, at minimum, the standard deviation of the 
#' variability that you desire.
#' You can also provide a mean value from which deviates will be created.
#' 
#' For \code{.harvest_periods} and {.harvest_offsets}, provide either a single integer
#' or a vector of a length that is a factor of \code{.n_patches}.
#'
#' @param .n_pops Number of aphid lines.
#' @param .n_patches Number of patches.
#' @param .harvest_periods Period of harvest cycle. Defaults to 
#'     \code{sappp::environ$cycle_length}.
#' @param .harvest_offsets Offset for harvest cycle. Defaults to 0.
#' @param aphid_density_0 Starting aphid densit(ies). Defaults to \code{NULL}.
#' @param aphid_surv_juv Aphid juvenile survival. Defaults to \code{NULL}.
#' @param aphid_surv_adult Vector of adult aphid survivals. Defaults to \code{NULL}.
#' @param aphid_repro Vector of aphid reproductive rates. Defaults to \code{NULL}.
#' @param K Aphid density dependency. Defaults to \code{NULL}.
#' @param disp_aphid Aphid dispersal rates. Defaults to \code{NULL}.
#' @param attack_surv Aphid wasp-attack survival rates, for singly and multiply 
#'     infected aphids. Defaults to \code{NULL}.
#' @param no_error Boolean for whether to include no process error. Defaults to 
#'     \code{TRUE}.
#' @param ... Other parameters associated with aphid lines and their infecting wasps.
#'     See \code{\link{all_pop_lists}} for info on how to provide these if you wish
#'     to vary by line or patch, and 
#'     see \code{vignette(parameters)} for the parameters available.
#'
#' @return
#' An SimPatches object.
#' 
#' @export
#' 
#' @examples
#' 
#' so <- make_sim_obj(
#'     .n_pops = 10, .n_patches = 8,
#'     .harvest_periods = 7,
#'     .harvest_offsets = c(0, 3),
#'     aphid_density_0 = NULL,
#'     aphid_surv_juv = list(.sd_pops = 2),
#'     aphid_surv_adult = list(.sd_pops = 1),
#'     aphid_repro = list(.sd_pops = 2),
#'     K = list(.sd_pops = 1),
#'     disp_aphid = list(.sd_pops = 1),
#'     attack_surv = NULL,
#'     no_error = TRUE)
#' sim <- so$simulate(630)
#' sim
#' 
#'
make_sim_obj <- function(
    .n_pops, .n_patches,
    .harvest_periods = NULL,
    .harvest_offsets = 0,
    aphid_density_0 = NULL,
    aphid_surv_juv = NULL,
    aphid_surv_adult = NULL,
    aphid_repro = NULL,
    K = NULL,
    disp_aphid = NULL,
    attack_surv = NULL,
    no_error = TRUE,
    ...) {
    
    if (is.null(.harvest_periods)) .harvest_periods <- sappp::environ$cycle_length
    
    # Get no-variance objects for those that are not input here
    # If they are input here, then get lists with variability
    par_programmed <- c('aphid_density_0', 'aphid_surv_juv', 'aphid_surv_adult', 
                        'aphid_repro', 'K', 'disp_aphid', 'attack_surv')
    for (.par in par_programmed) {
        if (is.null(eval(parse(text = .par)))) {
            # For .par = "K", the below paste0 call generates the following string:
            # "K <- do.call(make_n_values, list(\"K\", 1, 1, 0))"
            eval(parse(text = paste0(
                .par, ' <- ', 'do.call(make_n_values, list("', .par, '", 1, 1, 0))')))
        } else if (is.list(eval(parse(text = .par)))) {
            if (is.null(eval(parse(text = paste0(.par, '$.sd_pops'))))) {
                # "K <- do.call(make_n_values, c(.par_name = \"K\", .n_pops = .n_pops, K))"
                eval(parse(text = paste0(
                    .par, ' <- ', 'do.call(make_n_values, c(.par_name = "', .par, 
                    '", .n_pops = .n_pops, .n_patches = .n_patches, ', .par, '))')))
            }
        }
    }
    
    # Create lists specifying info for object creation
    parL <- all_pop_lists(.n_pops, .n_patches,
                          aphid_density_0 = aphid_density_0, 
                          aphid_surv_juv = aphid_surv_juv, 
                          aphid_surv_adult = aphid_surv_adult, 
                          aphid_repro = aphid_repro, 
                          K = K, 
                          disp_aphid = disp_aphid, 
                          attack_surv = attack_surv, ...)
    
    # Remove error if desired
    if (no_error) {
        for (i in 1:.n_patches) for (j in 1:.n_pops) {
            parL[[i]][[j]]$sigma_x <- 0
            parL[[i]][[j]]$sigma_y <- 0
            parL[[i]][[j]]$rho <- 0
            parL[[i]][[j]]$demog_mult <- 0
        }
    }
    
    # Creating vectors of harvest cycle lengths and offsets
    stopifnot((.n_patches %% length(.harvest_periods)) == 0)
    stopifnot((.n_patches %% length(.harvest_offsets)) == 0)
    .harvest_periods <- rep(.harvest_periods, .n_patches / length(.harvest_periods))
    .harvest_offsets <- rep(.harvest_offsets, .n_patches / length(.harvest_offsets))
    
    sp <- new(SimPatches, parL, .harvest_periods, .harvest_offsets)
    
    return(sp)
}




