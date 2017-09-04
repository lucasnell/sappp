library(Rcpp)
# devtools::clean_dll()
compileAttributes(); devtools::load_all()



# 
# # Estimated values from paper (symbol; value from paper)
# a <- 2.5        # parasitoid attack rate (a; 2.32)
# K <- 0.0005     # aphid density dependence (K; 0.000467)
# K_y <- 0.0006    # parasitized aphid density dependence (K_y; 0.00073319)
# k <- 0.1811    # aggregation parameter of the negative binomial distribution (k; 0.35)
# h <- 0.0363     # parasitoid attack rate handling time (h; 0.008, 0.029 at 27 deg C)
# s_y <- 0.55      # parasitoid adult daily survival (s_y; 0.69)
# sigma_x <- 0.44  # 0  # environmental standard deviation for aphids (sigma_x; 0.44)
# sigma_y <- 0.35  # 0  # environmental standard deviation for parasitoids (sigma_y; 0.70)
# rho <- 2 / (1 + exp(-sigma_y)) - 1  # environmental correlation among instars (rho; 1.0)
# 
# 
# sap::wasp_attack$a
# sap::populations$K
# sap::populations$K_y
# sap::wasp_attack$k
# sap::wasp_attack$h
# sap::populations$s_y




n_cycles <- 20
cycle_length <- 30
max_time <- cycle_length * (1 + n_cycles)
harvest_times <- list(c(cycle_length * 1:n_cycles),
                      c(cycle_length * 1, cycle_length * (2:(n_cycles-1)) - cycle_length/2,
                         cycle_length * n_cycles))
harvest_times <- lapply(harvest_times, function(x) head(x, -1))

n_pops = 2
n_patches = 2
pl <- all_pop_lists(n_pops, n_patches, 
                    a = 2.5, 
                    K = 0.0005, 
                    K_y = 0.0006, 
                    k = 0.1811, 
                    h = 0.0363, 
                    s_y = 0.55,
                    aphid_density_0 = list(
                        (1 - sap::populations$prop_resist) * sap::populations$aphids_0,
                        sap::populations$prop_resist * sap::populations$aphids_0),
                    # sigma_x = 0, sigma_y = 0, rho = 0, demog_mult = 0,
                    attack_surv = list("susceptible", "resistant"), 
                    aphid_surv_adult = list('high', 'low'),
                    aphid_repro = list('high', 'low'))
sp <- new(SimPatches, pl, harvest_times)
          # rep(sap::environ$cycle_length, n_patches), rep(c(0,15), n_pops / 2))
out <- sp$simulate(100, rng_seed = sample.int(2^31-1,1))


# Figure 1
{
    par(mar = c(b = 4, l = 2.5, t = 1, r = 1))
    # X is sum of all living aphids
    Xr <- out$aphids[,2,] + out$parasit[,2,]
    Xs <- out$aphids[,1,] + out$parasit[,1,]
    # Y is proportion of live aphids that are parasitized
    Yr <- out$parasit[,2,] / Xr
    Ys <- out$parasit[,1,] / Xs
    Ymax <- 1 * max(c(Xr[,1]+Xs[,1], Xr[,2]+Xs[,2]))
    min_time <- 1
    max_time <- min(630, nrow(Xr))
    times <- 1:(max_time-min_time+1)
    plot(times,Xr[min_time:max_time,1]+Xs[min_time:max_time,1],type = 'l', 
         col = 'dodgerblue', ylab = '', xlab = 'time', main = '', ylim = c(1,Ymax))
    lines(times,Xr[min_time:max_time,2]+Xs[min_time:max_time,2], 
          col = 'dodgerblue', lty = 2)
    lines(times,Ymax * (Yr[min_time:max_time,1]+Ys[min_time:max_time,1])/2, 
          col = 'firebrick')
    lines(times,Ymax * (Yr[min_time:max_time,2]+Ys[min_time:max_time,2])/2, 
          col = 'firebrick', lty = 2)
    lines(times,Ymax*rowSums(Xr[min_time:max_time,])/(
        rowSums(Xr[min_time:max_time,]) + rowSums(Xs[min_time:max_time,])), 
        col = 'black')
}
# Blue is aphid abundances
# Red is parasitoid abundances
# Black is (resistant aphids) / (total aphids)
# Dotted lines are field # 2 (different harvesting regime)



#



# AphidWasp_one_patch()








# -------
# Susceptible aphid line info
# -------
sus_line <- make_AphidWasp(
    "susceptible",
    aphid_density_0 = (1 - sap::populations$prop_resist) * sap::populations$aphids_0)
sus_line



# -------
# Resistant aphid line info
# -------
# Difference is resistance, lower reproduction, and lower starting density
res_line <- make_AphidWasp(
    "resistant",
    attack_surv = sap::wasp_attack$attack_surv,
    aphid_density_0 = sap::populations$prop_resist * sap::populations$aphids_0,
    aphid_repro = sap::populations$repro$low)
res_line


sourceCpp(code = 
"
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//[[Rcpp::export]]
void test(uint rows, uint cols, uint slices, uint t, arma::mat m) {
    arma::cube x(rows, cols, slices, arma::fill::randu);
    x.print(Rcout);
    // x.slice(s).print(Rcout);
    for (arma::uword i = 0; i < x.n_slices; i++) x.slice(i).row(t) = m.row(i);
    x..tube(arma::span(t, t), arma::span()).print(Rcout);
    return;
}
")

test(4, 3, 3, 1, matrix(11:19, 3, 3, byrow = TRUE))


array(0, c(4, 3, 2))


which(sapply(1:100, do_harvest))






# =================================================================================
# =================================================================================

#           APHID-LINE INFO

# =================================================================================
# =================================================================================

# Constant info about an aphid line (and wasps that parasitize them)
# This should be used among fields for the same aphid line
aphid_const <- setRefClass(
    
    'aphid_const', 
    
    fields = list(
        leslie = 'matrix',          # Leslie matrix with survival and reproduction
        X_0 = 'numeric',            # initial aphid abundances by stage
        Y_0 = 'numeric',            # initial wasp abundances by stage
        rel_attack = 'numeric',     # relative wasp attack rates by aphid stage
        sex_ratio = 'numeric',      # proportion female
        a = 'numeric',              # overall parasitoid attack rate
        K = 'numeric',              # aphid density dependence
        K_y = 'numeric',            # parasitized aphid density dependence
        k = 'numeric',              # aggregation parameter of the nbinom distribution
        h = 'numeric',              # parasitoid attack rate handling time
        s_y = 'numeric',            # parasitoid adult daily survival
        sigma_x = 'numeric',        # environmental standard deviation for aphids
        sigma_y = 'numeric',        # environmental standard deviation for parasitoids
        rho = 'numeric',            # environmental correlation among instars
        attack_surv = 'numeric',    # survival rates of singly & multiply attacked aphids
        harvest_surv = 'numeric',   # survival rate for aphids during a harvest
        disp_aphid = 'numeric',     # dispersal rate for aphids
        disp_stages = 'integer',    # stages in which dispersal occurs in aphids
        disp_wasp = 'numeric',      # dispersal rate for wasps
        pred_rate = 'numeric',      # predation on aphids and mummies
        n_aphid_stages = 'integer', # number of aphid stages (i.e., days)
        n_wasp_stages = 'integer',  # number of wasp stages (i.e., days)
        mum_days = 'integer'        # number of days per mummy stage (aphid alive & dead)
    )

)

aphid_const$methods(
    
    initialize = function(leslie, rel_attack, 
                          a, K, K_y, k, h, s_y, sigma_x, sigma_y, rho,
                          disp_stages, mum_days, 
                          aphid_density_0, wasp_density_0 = 1, 
                          sex_ratio = 0.5, attack_surv = rep(0.0,2), harvest_surv = 0.05,
                          disp_aphid = 0.05, disp_wasp = 1,
                          pred_rate = 0.8) {
        'Method for initializing an object of aphid_const class'
        # --------
        # Checking for nonsense
        # --------
        stopifnot(length(a) == 1, length(K) == 1, length(K_y) == 1, length(k) == 1, 
                  length(h) == 1, length(s_y) == 1, length(sigma_x) == 1,  
                  length(sigma_y) == 1,  length(rho) == 1, length(wasp_density_0) == 1, 
                  length(aphid_density_0) == 1, length(sex_ratio) == 1, 
                  length(harvest_surv) == 1, length(disp_aphid) == 1, 
                  length(disp_wasp) == 1, length(pred_rate) == 1)
        stopifnot(length(mum_days) == 2)
        
        # --------
        # Fields directly provided to function
        # --------
        leslie <<- leslie
        rel_attack <<- as.numeric(rel_attack)
        a <<- as.numeric(a)
        K <<- as.numeric(K)
        K_y <<- as.numeric(K_y)
        k <<- as.numeric(k)
        h <<- as.numeric(h)
        s_y <<- as.numeric(s_y)
        sigma_x <<- as.numeric(sigma_x)
        sigma_y <<- as.numeric(sigma_y)
        rho <<- as.numeric(rho)
        if (length(attack_surv) >= 2) {
            attack_surv <<- as.numeric(attack_surv)
        } else {
            attack_surv <<- numeric(2)
        }
        disp_stages <<- as.integer(disp_stages)
        mum_days <<- as.integer(mum_days)
        sex_ratio <<- as.numeric(sex_ratio)
        harvest_surv <<- as.numeric(harvest_surv)
        disp_aphid <<- as.numeric(disp_aphid)
        disp_wasp <<- as.numeric(disp_wasp)
        pred_rate <<- as.numeric(pred_rate)
        
        # --------
        # Fields calculated from inputs
        # --------
        n_aphid_stages <<- as.integer(nrow(leslie))
        n_wasp_stages <<- as.integer(sum(mum_days) + 1)
        Y_0 <<- as.numeric(c(rep(0, sum(mum_days)), wasp_density_0))
        X_0 <<- as.numeric(aphid_density_0 * leslie_sad(leslie))

    },
    
    show = function() {
        'Method for showing an object of aphid_const class'
        cat("<< Aphid line constant info >> \n")
        cat("Parasitoid resistance vector: ", attack_surv, "\n")
        cat("First 7x7 cells of the Leslie matrix:\n")
        print(leslie[1:7,1:7])
    }
    
)




# Info about an aphid line (and wasps that parasitize them), some of which changes 
# through time.
# This should be used for a single line in a single field.
# The line_info field can use shared memory for all of the same line across fields
aphid_iter <- setRefClass(
    
    'aphid_iter',
    
    fields = list(
        line_info = 'ANY',          # Aphid line info (doesn't change through time)
        X_t = 'numeric',            # Aphid density at time t
        X_t1 = 'numeric',           # Aphid density at time t+1
        Y_t = 'numeric',            # Wasp density at time t
        Y_t1 = 'numeric',           # Wasp density at time t+1
        A = 'numeric'               # Attack probabilities at time t+1
    )
)

aphid_iter$methods(
    
    initialize = function(line_info) {
        'Method for initializing an object of aphid_iter class'
        stopifnot('aphid_const' %in% class(line_info))
        
        line_info <<- line_info
        
        X_t <<- line_info$X_0
        X_t1 <<- line_info$X_0
        Y_t <<- line_info$Y_0
        Y_t1 <<- line_info$Y_0
        A <<- numeric(length(X_t))
    },
    
    show = function() {
        'Method for showing an object of aphid_const class'
        cat("<< Aphid line changing info >> \n")
        cat("Parasitoid resistance vector: ", line_info$attack_surv, "\n")
        cat("First 7x7 cells of the Leslie matrix:\n")
        print(line_info$leslie[1:7,1:7])
    }
    
)





# =================================================================================
# =================================================================================

#           PATCH INFO

# =================================================================================
# =================================================================================

# Info about a patch, some of which changes through time
patch <- setRefClass(
    
    'patch', 
    
    fields = list(
        harvest_times = 'integer',  # times when harvesting occurs
        max_time = 'integer',       # number of time points to simulate
        z = 'numeric',              # Sum of all living aphids at time t
        x = 'numeric',              # Sum of non-parasitized aphids at time t
        Y_m =  'numeric'            # Total number of adult wasps
    )
)


patch$methods(
    
    initialize = function(harvest_times, max_time) {
        'Method for initializing an object of patch class'
        stopifnot(length(max_time) == 1)
        
        harvest_times <<- harvest_times
        max_time <<- max_time
        z <<- numeric(1)
        x <<- numeric(1)
        Y_m <<- numeric(1)
    },
    
    # Density dependence for aphids (note that equation is different in paper)
    S = function(K) {
        1 / (1 + K * z)
    },
    
    # Density dependence for wasps (note that equation is different in paper)
    S_y = function(K_y) {
        1 / (1 + K_y * z)
    },
    
    
    # Left off --> not sure the below function works!!
    
    # Run at the beginning of each iteration to have updated values
    update = function(aphid_obj_list) {
        tmp_list <- lapply(aphid_obj_list, 
                           function(a) {
                               linf <- a$line_info
                               cbind(sum(a$X_t[,]), sum(a$Y_t[1:(linf$mum_days[1])]),
                                     sum(a$Y_t[length(a$Y_t)]))
                           })
        sum_mat <- do.call(rbind, tmp_list)
        z <<- sum(sum_mat[,1:2])
        x <<- sum(sum_mat[,1])
        Y_m <<- sum(sum_mat[,3])
    },
    
    show = function() {
        'Method for showing an object of aphid_const class'
        cat("<< Aphid line changing info >> \n")
        cat("Parasitoid resistance vector: ", line_info$attack_surv, "\n")
        cat("First 7x7 cells of the Leslie matrix:\n")
        print(line_info$leslie[1:7,1:7])
    }
    
)







# =================================================================================
# =================================================================================

#           FULL SIMULATION INFO

# =================================================================================
# =================================================================================


# Combining the classes above to make one object per simulation
