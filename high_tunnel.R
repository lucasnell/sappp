# 
# This script seeks to replicate the high-tunnel model that Tony sent me Aug 2017
# 
# by:  Lucas A. Nell
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix) # for sparse matrices

# Provides functions instar_to_stage, leslie_matrix, leslie_sad, attack_probs, 
# parasitoid_abunds
Rcpp::sourceCpp('high_tunnel.cpp')


base_plot <- function() {
    Ymax <- 1.2 * max(Xr+Xs)
    Ymax2 <- 1.2 * max(Xr+Xs)
    min_time <- 1
    
    # Figure 1
    plot(1:(max_time-min_time+1),Xr[min_time:max_time,1]+Xs[min_time:max_time,1],
         type = 'l', col = 'dodgerblue', ylab = '', xlab = 'time', main = '')
    lines(1:(max_time-min_time+1),Xr[min_time:max_time,2]+Xs[min_time:max_time,2],
          col = 'dodgerblue', lty = 2)
    lines(1:(max_time-min_time+1),Ymax2 * (Yr[min_time:max_time,1]+Ys[min_time:max_time,1]),
          col = 'firebrick')
    lines(1:(max_time-min_time+1),Ymax2 * (Yr[min_time:max_time,2]+Ys[min_time:max_time,2]),
          col = 'firebrick', lty = 2)
    lines(1:(max_time-min_time+1),Ymax2*rowSums(Xr[min_time:max_time,])/
              (rowSums(Xr[min_time:max_time,]) + rowSums(Xs[min_time:max_time,])),
          col = 'black')
    
    # Blue is aphid abundances
    # Red is parasitoid abundances
    # Black is (resistant aphids) / (total aphids)
    # Dotted lines are field # 2 (different harvesting regime)
}



# process error for aphids, parasitized aphids, and adult parasitoids

# mat X_t1 = xtr
# mat Y_t1 = ytr
# mat X_t = xr[,i]
# mat Y_t = yr[,i]
# double sigma_x = sigma_x
# double sigma_y = sigma_y
# double rho = rho
# double z = z
# double Y_m = yr[nrow(yr),i]
# uword total_stages = n_aphid_stages+n_wasp_stages
#   Total days for living aphids (i.e., not parasitized or parasitized 
#   but not yet a mummy):
# uword living_aphids = n_aphid_stages + mum_days[1]
# double sigma_d_mult = 1  # A way to scale the demographic var. (for wasps and aphids)


# mat X_t1, mat Y_t1, mat X_t, mat Y_t, double sigma_x, double sigma_y,
# double rho, double z, double Y_m, uword total_stages, uword living_aphid



process_error <- function(X_t1, Y_t1, X_t, Y_t, sigma_x, sigma_y, rho, z, Y_m, 
                          total_stages, living_aphids, sigma_d_mult = 1) {
    
    Se = matrix(0, total_stages, total_stages);
    
    if (sigma_d_mult == 0 & sigma_x == 0 & sigma_y == 0) {
        return(list(aphids = X_t1, wasps = Y_t1));
    }
    
    # Aphid (both parasitized and not) process error
    Se[1:living_aphids,1:living_aphids] =
        # Version from paper:
        (sigma_x^2 + sigma_d_mult * min(0.5, 1 / abs(1 + z))) *
        (rho * matrix(1,living_aphids,living_aphids) + (1-rho) * diag(1,living_aphids));
        # # Version Tony sent:
        # sigma_x^2 * (rho * matrix(1,living_aphids,living_aphids) + (1-rho) * 
        # diag(1,living_aphids))
    
    # Mummy process error, turning back to zero
    Se[(living_aphids+1):(nrow(Se)-1),(living_aphids+1):(ncol(Se)-1)] = 0;
    Se[1:living_aphids,(living_aphids+1):(ncol(Se)-1)] = 0;
    Se[(living_aphids+1):(nrow(Se)-1),1:living_aphids] = 0;
    
    # Adult parasitoid process error
    Se[nrow(Se),ncol(Se)] = sigma_y^2 + sigma_d_mult * min(0.5, 1 / abs(1 + Y_m));
    
    # chol doesn't work with zeros on diagonal
    non_zero = which(diag(Se) > 0);
    
    # Cholesky decomposition of Se so output has correct variance-covariance matrix
    #   "a vector of independent normal random variables,
    #   when multiplied by the transpose of the Cholesky deposition of [Se] will
    #   have covariance matrix equal to [Se]."
    chol_decomp = t(chol(Se[non_zero,non_zero]));
    
    # Random numbers from distribution N(0,1)
    rnd = rnorm(length(non_zero));
    
    # Making each element of rnd have correct variance-covariance matrix
    E = chol_decomp %*% rnd;
    
    nz_aphid = non_zero[non_zero <= nrow(X_t1)];
    nz_wasp = non_zero[non_zero > nrow(X_t1)] - nrow(X_t1);
    
    X_t1_e = X_t1;
    Y_t1_e = Y_t1;
    
    X_t1_e[nz_aphid,] = X_t1_e[nz_aphid,] * exp(E[1:length(nz_aphid)]);
    Y_t1_e[nz_wasp,] = Y_t1_e[nz_wasp,] * exp(E[(length(nz_aphid)+1):nrow(E)]);
    
    # Because we used normal distributions to approximate demographic and environmental 
    # stochasticity, it is possible for aphids and parasitoids to 
    # "spontaneously appear" when the estimate of e(t) is large. To disallow this 
    # possibility, the number of aphids and parasitized aphids in a given age class 
    # on day t was not allowed to exceed the number in the preceding age class on 
    # day t – 1.
    
    for (i in 2:length(X_t1_e)) {
        X_t1_e[i] = min(X_t1_e[i], X_t[(i-1)]);
        # if (X_t1_e[i] > X_t[(i-1)]) {
        #     X_t1_e[i] = X_t[(i-1)];
        # }
    }
    # Not going to the end for parasitoids bc you can have more adults than mummies
    # bc adults stay in that stage for multiple days
    for (i in 2:(length(Y_t1_e)-1)) {
        Y_t1_e[i] = min(Y_t1_e[i], Y_t[(i-1)]);
        # if (Y_t1_e[i] > Y_t[(i-1)]) {
        #     Y_t1_e[i] = Y_t[(i-1)];
        # }
    }
    
    return(list(aphids = X_t1_e, wasps = Y_t1_e));
}


# NOTE: The code is set up for 2 clones that have different life histories.
# The life histories are taken from lab experiments and correspond to
# demographic parameters at 20 and 27 C. See Meisner et al. 2014.
# 
# The juvenile and adult demographic rates are separated, so you can pick
# them separately for the different clones. "Resistant" clones are
# resistant to parasitism.
# 
# rows: 1=resistant, 2=susceptible
# columns: 1=aphid juv, 2=aphid adult
# values 1=low growth rate, 2=high growth rate
clone <- rbind(c(2, 1), c(2, 2))
# Values dictate which rows will be selected from matrices of demographic rates below
# (stage_days, surv_juv, surv_adult, repro)



# =============================================
# lab parameters
# =============================================

n_aphid_stages <- 32
n_lines <- 2


# Number of days per instar
stage_days <- rbind(c(2, 2, 2, 2, 19), c(1, 1, 1, 2, 23))
# Number of days for parasitized (but still living) aphids and mummies
mum_days <- cbind(7, 3)

n_wasp_stages <- sum(mum_days) + 1

# juvenile survival
surv_juv <- rbind(0.9745, 0.9849)

# adult survival
surv_adult <- rbind(
    c(1.0000, 0.9949, 0.9818, 0.9534, 0.8805, 0.8367, 0.8532, 0.8786, 0.8823, 
      0.8748, 0.8636, 0.8394, 0.8118, 0.8096, 0.8240, 0.8333, 0.7544, 0.5859, 
      0.4155, 0.2216),
    c(1.0000, 0.9986, 0.9951, 0.9874, 0.9675, 0.9552, 0.9550, 0.9549, 0.9462, 
      0.8992, 0.8571, 0.8408, 0.8281, 0.8062, 0.7699, 0.7500, 0.7559, 0.7649, 
      0.7240, 0.4367))
surv_adult <- cbind(surv_adult, matrix(0,n_lines,180))

# reproduction
repro <- rbind(
    c(0, 2.5925, 4.4312, 5.1403, 5.5190, 5.6633, 5.6010, 5.4577, 5.2904, 5.0613, 4.6970, 
      3.3577, 1.5946, 1.0817, 0.9666, 0.8333, 0.4689, 0.0709, 
      0, 0, 0, 0),
    c(0, 3.1975, 5.4563, 6.2996, 6.7372, 6.9030, 6.8210, 6.6100, 6.1962, 5.1653, 4.1837, 
      3.6029, 3.1023, 2.4799, 1.6909, 1.1750, 1.0148, 0.9096, 
      0.7821, 0.6430, 0.5000, 0.3531)
)
repro <- cbind(repro, matrix(0,n_lines,178))



# Relative attack rate on the different instars from Ives et al 1999
# `instar_to_stage` converts these values from per-instar to per-day
rel_attack <- instar_to_stage(cbind(0.12, 0.27, 0.39, 0.16, 0.06), 
                              n_aphid_stages, stage_days)




sex_ratio <- 0.5

# Estimated values from paper (symbol; value from paper)
a <- 2.5        # parasitoid attack rate (a; 2.32)
K <- 0.0005     # aphid density dependence (K; 0.000467)
K_y <- 0.0006    # parasitized aphid density dependence (K_y; 0.00073319)
k <- 0.1811    # aggregation parameter of the negative binomial distribution (k; 0.35)
h <- 0.0363     # parasitoid attack rate handling time (h; 0.008, 0.029 at 27 deg C)
s_y <- 0.55      # parasitoid adult daily survival (s_y; 0.69)
sigma_x <- 0.44  # 0  # environmental standard deviation for aphids (sigma_x; 0.44)
sigma_y <- 0.35  # 0  # environmental standard deviation for parasitoids (sigma_y; 0.70)
rho <- 2 / (1 + exp(-sigma_y)) - 1  # environmental correlation among instars (rho; 1.0)




# not used at all:
# mum_detect <- 0.3923

# For measurement error part:
# sm <- 33.5455

# These are the survivals of singly attacked and multiply attacked
# resistant aphids
resist_surv <- cbind(0.9, 0.6)




# =============================================
# set up Leslie matrices
# =============================================

# resistant clones
# ------
leslie_r <- leslie_matrix(n_aphid_stages, stage_days, clone[1,], surv_juv, 
                          surv_adult, repro)
sad_r <- leslie_sad(leslie_r)


# susceptible clones
# ------
leslie_s <- leslie_matrix(n_aphid_stages, stage_days, clone[2,], surv_juv, 
                          surv_adult, repro)
sad_s <- leslie_sad(leslie_s)



# =============================================
# field parameters: This is set up to have different harvesting patterns
# between n_fields fields.
# =============================================

# number of fields
n_fields <- 2

# kill rate at harvesting
kill <- 0.05

# dispersal rates between fields for aphids, adult wasps
disp_aphid <- 0.05
disp_wasp <- 1
# Predation rate
pred_rate <- 0.8

# initial densities of aphids and parasitoids
init_x <- 20
init_y <- 1

# time between harvests
n_cycles <- 20
cycle_length <- 30

# Proportion of resistant clones
prop_resist <- 0.05



# =============================================
# run program
# =============================================

# Initial densities of aphids by stage
xr <- prop_resist * init_x * sad_r %*% matrix(1,1,n_fields)
xs <- (1-prop_resist) * init_x * sad_s %*% matrix(1,1,n_fields)

# Initial parasitoid densities by stage (starting with no parasitized aphids or mummies)
yr <- init_y * rbind(matrix(0, sum(mum_days), n_fields), c(1, 1))
ys <- init_y * rbind(matrix(0, sum(mum_days), n_fields), c(1, 1))

# Setting total time (days) and times for harvesting
max_time <- cycle_length * (1 + n_cycles)
harvest_times <- rbind(c(cycle_length * 1:n_cycles),
                      c(cycle_length * 1, cycle_length * (2:(n_cycles-1)) - cycle_length/2,
                        cycle_length * n_cycles))














HighTunnelExptSimfunct <- function(
    xr,xs,yr,ys,a,resist_surv,K,K_y,k,h,s_y,
    sigma_x,sigma_y,rho,sigma_d_mult,
    # sm,
    max_time,n_fields,kill,
    harvest_times,disp_aphid,disp_wasp,pred_rate,
    n_aphid_stages, n_wasp_stages, stage_days, mum_days, surv_juv, surv_adult,
    rel_attack, sex_ratio, leslie_r, leslie_s, clone) {
    
    # Storing aphid (X) and wasp (Y) for resistant (r) and susceptible (s) clones in 
    # n_fields fields
    Xr <- matrix(0, max_time, n_fields)
    Yr <- matrix(0, max_time, n_fields)
    Xs <- matrix(0, max_time, n_fields)
    Ys <- matrix(0, max_time, n_fields)
    
    # Leslie matrices for resistant and susceptible aphids
    LLr <- Matrix(leslie_r, sparse = TRUE)
    LLs <- Matrix(leslie_s, sparse = TRUE)
    
    for (t in 1:max_time) {
        for (i in 1:n_fields) {
    # t=1;i=1
            
            # Total number of living aphids (z in the paper)
            z <- sum(c(xr[,i], xs[,i], yr[1:mum_days[1], i], ys[1:mum_days[1], i]))
            
            # Equivalent to S(z(t)) and S_y(z(t)) [no idea why equations are different]
            St <- 1 / (1 + K * z)
            S_yt <- 1 / (1 + K_y * z)
            
            # Matrices of attack probabilities using equation 6 from paper
            As <- attack_probs(a = a, p_i = rel_attack[clone[1,1],], 
                               Y_m = yr[nrow(yr),i] + ys[nrow(ys),i], 
                               x = z, h = h, k = k, resist_surv = numeric(0))
            Ar <- attack_probs(a = a, p_i = rel_attack[clone[2,1],], 
                               Y_m = yr[nrow(yr),i] + ys[nrow(ys),i], 
                               x = z, h = h, k = k, resist_surv = resist_surv)
            
            # X(t+1) (first line of equation 2; pred_rate was added later)
            xtr <- (pred_rate * St * Ar) * as.matrix(LLr %*% xr[,i])
            xts <- (pred_rate * St * As) * as.matrix(LLs %*% xs[,i])
            
            # Filling in column of parasitoid stage abundances
            ytr <- parasitoid_abunds(S_y_zt = S_yt, A = Ar, L = LLr, X = xr[,i], 
                                     Y_t = yr[,i], s_i = surv_juv[clone[1,1]], 
                                     s_y = s_y, m_1 = mum_days[1], sex_ratio = sex_ratio, 
                                     pred_rate = pred_rate)
            yts <- parasitoid_abunds(S_y_zt = S_yt, A = As, L = LLs, X = xs[,i], 
                                     Y_t = ys[,i], s_i = surv_juv[clone[2,1]], 
                                     s_y = s_y, m_1 = mum_days[1], sex_ratio = sex_ratio, 
                                     pred_rate = pred_rate)
            
            
            # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # # This code for stochasticity is dead
            # # process error for aphids, parasitized aphids, and adult parasitoids
            # Se <- matrix(0, n_aphid_stages+n_wasp_stages, n_aphid_stages+n_wasp_stages)
            # 
            # # Total days for living aphids (i.e., not parasitized or parasitized 
            # # but not yet a mummy)
            # living_aphids <- n_aphid_stages + mum_days[1]
            # 
            # # Aphid (both parasitized and not) process error
            # Se[1:living_aphids,1:living_aphids] <-
            #     # Version from paper:
            #     (sigma_x^2 + min(0.5, 1 / abs(1 + z))) *
            #     (rho * matrix(1,living_aphids,living_aphids)+(1-rho) * diag(1,living_aphids));
            #     # # Version Tony sent:
            # 	# sigma_x^2 * (rho * matrix(1,living_aphids,living_aphids) + (1-rho) * diag(1,living_aphids))
            # 
            # # Mummy process error (this isn't really needed, but I like it here for 
            # # transparency)
            # Se[(living_aphids+1):(nrow(Se)-1),(living_aphids+1):(ncol(Se)-1)] <- 0
            # Se[1:living_aphids,(living_aphids+1):(ncol(Se)-1)] <- 0
            # Se[(living_aphids+1):(nrow(Se)-1),1:living_aphids] <- 0
            # 
            # # Adult parasitoid process error
            # Se[nrow(Se),ncol(Se)] <- sigma_y^2 + min(0.5, 1 / abs(1 + yr[nrow(yr),i]));
            # ypick <- which(diag(Se) > 0)
            # 
            # iD <- t(chol(Se[ypick,ypick]))
            # E <- iD %*% rnorm(length(ypick))
            # 
            # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            error_list <- process_error(X_t1 = xtr, Y_t1 = ytr,
                                        X_t = xr[,i], Y_t = yr[,i],
                                        sigma_x = sigma_x, sigma_y = sigma_y, rho = rho,
                                        z = z, Y_m = yr[nrow(yr),i],
                                        total_stages = n_aphid_stages + n_wasp_stages,
                                        living_aphids = n_aphid_stages + mum_days[1],
                                        sigma_d_mult = sigma_d_mult)

            xtr <- error_list$aphids
            ytr <- error_list$wasps

            error_list <- process_error(X_t1 = xts, Y_t1 = yts,
                                        X_t = xs[,i], Y_t = ys[,i],
                                        sigma_x = sigma_x, sigma_y = sigma_y, rho = rho,
                                        z = z, Y_m = ys[nrow(ys),i],
                                        total_stages = n_aphid_stages + n_wasp_stages,
                                        living_aphids = n_aphid_stages + mum_days[1],
                                        sigma_d_mult = sigma_d_mult)

            xts <- error_list$aphids
            yts <- error_list$wasps
            
            if (t %in% harvest_times[i,1:(ncol(harvest_times)-1)]) {
                # Kill non-parasitized aphids
                xtr <- kill * xtr
                xts <- kill * xts
                # Kill parasitized (but still living) aphids at the same rate
                ytr[1:mum_days[1]] <- kill * ytr[1:mum_days[1]]
                yts[1:mum_days[1]] <- kill * yts[1:mum_days[1]]
                # Kill all mummies
                ytr[(mum_days[1]+1):(length(ytr)-1)] <- 0
                yts[(mum_days[1]+1):(length(yts)-1)] <- 0
            }
            
            # Filling in values for the next iteration (xr, xs, yr, ys hold time t info)
            xr[,i] <- xtr
            xs[,i] <- xts
            yr[,i] <- ytr
            ys[,i] <- yts
        }
        
        # First 4 instars apparently don't disperse
        dispersing <- disp_aphid * cbind(rowMeans(xr[(sum(stage_days[clone[1,1],1:4])+1):nrow(xr),]))
        xr[(sum(stage_days[clone[1,1],1:4])+1):nrow(xr),] <- (1-disp_aphid) *
            xr[(sum(stage_days[clone[1,1],1:4])+1):nrow(xr),] + 
            dispersing %*% matrix(1,1,n_fields)
        
        dispersing <- disp_aphid * cbind(rowMeans(xs[(sum(stage_days[clone[2,1],1:4])+1):nrow(xs),]))
        xs[(sum(stage_days[clone[2,1],1:4])+1):nrow(xs),] <- (1-disp_aphid) *
            xs[(sum(stage_days[clone[2,1],1:4])+1):nrow(xs),] + 
            dispersing %*% matrix(1,1,n_fields)
        
        dispersingw <- disp_wasp * mean(yr[nrow(yr),])
        yr[nrow(yr),] <- ((1-disp_wasp) * yr[nrow(yr),] + dispersingw)
        
        dispersingw <- disp_wasp * mean(ys[nrow(ys),])
        ys[nrow(ys),] <- ((1-disp_wasp) * ys[nrow(ys),] + dispersingw)
        
        # X is sum of all living aphids
        Xr[t,] <- colSums(xr) + colSums(yr[1:mum_days[1],])
        Xs[t,] <- colSums(xs) + colSums(ys[1:mum_days[1],])
        
        # Y is proportion of live aphids that are parasitized
        Yr[t,] <- colSums(yr[1:mum_days[1],]) / Xr[t,]
        Ys[t,] <- colSums(ys[1:mum_days[1],]) / Xs[t,]
        
    }
    

    out_list <- list(Xr = Xr, Xs = Xs, Yr = Yr, Ys = Ys)

    return(out_list)
}



sigma_d_mult = 1

set.seed(60253704)
out_list <- HighTunnelExptSimfunct(xr,xs,yr,ys,a,resist_surv,K,K_y,k,h,s_y,
                                   sigma_x,sigma_y,rho,sigma_d_mult, # <- Full error
                                   # 0,0,0,1,  # <- No environmental error
                                   # sigma_x,sigma_y,rho,0, # <- No demographic error
                                   # 0,0,0,0, # <- No error
                                   max_time,n_fields,kill,harvest_times,disp_aphid,
                                   disp_wasp,pred_rate,
                                   n_aphid_stages, n_wasp_stages, stage_days, mum_days, 
                                   surv_juv, surv_adult, 
                                   rel_attack, sex_ratio, leslie_r, leslie_s, clone)
# Assigning out_list values to global ones for objects Xr, Xs, Yr, Ys
Xr <- out_list$Xr
Xs <- out_list$Xs
Yr <- out_list$Yr
Ys <- out_list$Ys

base_plot()




aphids <- as_data_frame(cbind(Xr, Xs)) %>%
    mutate(time = 1:nrow(Xs), `1` = V1+V3, `2` = V2+V4,
           r_prop = (V1+V2)/(V1+V2+V3+V4)) %>% # <-- (resistant aphids) / (total aphids)
    select(-1:-4) %>% 
    gather('field', 'density', 2:3, factor_key = TRUE)

wasps <- as_data_frame(cbind(Ys, Yr)) %>%
    mutate(time = 1:nrow(Ys), `1` = (V1+V3)/2, `2` = (V2+V4)/2) %>% 
    select(-1:-4) %>% 
    gather('field', 'density', 2:3, factor_key = TRUE)

# axis_mult = max(aphids$density)
axis_mult = 330

ggplot(aphids, aes(time)) +
    theme_classic() +
    theme(legend.position = c(0.01, 1), legend.direction = 'horizontal',
          legend.justification = c('left', 'top'), legend.margin = margin(0,0,0,0),
          legend.background = element_rect(fill = NA)) +
    geom_line(aes(y = density, linetype = field), color = 'dodgerblue') +
    geom_line(data = wasps, aes(y = density * axis_mult, linetype = field), 
              color = 'firebrick') +
    geom_line(aes(y = r_prop * axis_mult)) +
    geom_text(data = data_frame(time = c(200, 350, 500), 
                                y = rep(axis_mult * 1.1, 3), 
                                lab = c('aphids', '% parasit.', '% resist.')),
              aes(y = y, label = lab), color = c('dodgerblue', 'firebrick', 'black'),
              hjust = c(0, 0.5, 1), vjust = 1) +
    scale_y_continuous('aphid density', limits = c(0, axis_mult * 1.1),
                       sec.axis = sec_axis(~ . * 100 / axis_mult, name = '%'))



base_plot()
