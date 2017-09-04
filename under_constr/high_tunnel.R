# 
# This script seeks to replicate the high-tunnel model that Tony sent me Aug 2017
# 
# by:  Lucas A. Nell
# 


library(dplyr)
library(tidyr)
library(ggplot2)
# library(Matrix) # for sparse matrices

# Provides functions instar_to_stage, leslie_matrix, leslie_sad, attack_probs, 
# parasitoid_abunds
Rcpp::sourceCpp('high_tunnel.cpp')

source('parameters.R')









HighTunnelExptSimfunct <- function(
    X_0r,X_0s,Y_0r,Y_0s,a,attack_surv,K,K_y,k,h,s_y,
    sigma_x,sigma_y,rho,sigma_d_mult,
    # sm,
    max_time,n_fields,harvest_surv,
    harvest_times,disp_aphid,disp_wasp,pred_rate,
    n_aphid_stages, n_wasp_stages, instar_days, mum_days, surv_juv, surv_adult,
    rel_attack, sex_ratio, leslie_r, leslie_s, clone) {
    
    # Storing aphid (X) and wasp (Y) for resistant (r) and susceptible (s) clones in 
    # n_fields fields
    Xr <- matrix(0, max_time, n_fields)
    Yr <- matrix(0, max_time, n_fields)
    Xs <- matrix(0, max_time, n_fields)
    Ys <- matrix(0, max_time, n_fields)
    
    # Leslie matrices for resistant and susceptible aphids
    LLr <- leslie_r
    LLs <- leslie_s
    
    for (t in 1:max_time) {
        for (i in 1:n_fields) {

            # Total number of living aphids (z in the paper)
            z <- sum(c(X_0r[,i], X_0s[,i], Y_0r[1:mum_days[1], i], Y_0s[1:mum_days[1], i]))
            
            # Total number of non-parasitized aphids (x in paper)
            x <- sum(c(X_0r[,i], X_0s[,i]))
            
            # Equivalent to S(z(t)) and S_y(z(t)) [no idea why equations are different]
            St <- 1 / (1 + K * z)
            S_yt <- 1 / (1 + K_y * z)
            
            # Total number of adult wasps
            Y_m <- Y_0r[nrow(Y_0r),i] + Y_0s[nrow(Y_0s),i]
            
            # Matrices of attack probabilities using equation 6 from paper
            # iterate_A(double Y_m, double x)
            Ar <- attack_probs(a = a, p_i = rel_attack[[clone[1,1]]], 
                               Y_m = Y_m, x = x, h = h, k = k, 
                               attack_surv = attack_surv)
            As <- attack_probs(a = a, p_i = rel_attack[[clone[2,1]]], 
                               Y_m = Y_m, x = x, h = h, k = k, 
                               attack_surv = numeric(0))
            
            # X(t+1) (first line of equation 2; pred_rate was added after paper)
            # iterate_X(double S)
            xtr <- (pred_rate * St * Ar) * as.matrix(LLr %*% X_0r[,i])
            xts <- (pred_rate * St * As) * as.matrix(LLs %*% X_0s[,i])
            
            # Filling in column of parasitoid stage abundances
            # iterate_Y(double S_y)
            ytr <- parasitoid_abunds(S_y_zt = S_yt, A = Ar, L = LLr, X = X_0r[,i], 
                                     Y_t = Y_0r[,i], s_i = surv_juv[[clone[1,1]]], 
                                     s_y = s_y, m_1 = mum_days[1], sex_ratio = sex_ratio, 
                                     pred_rate = pred_rate)
            yts <- parasitoid_abunds(S_y_zt = S_yt, A = As, L = LLs, X = X_0s[,i], 
                                     Y_t = Y_0s[,i], s_i = surv_juv[[clone[2,1]]], 
                                     s_y = s_y, m_1 = mum_days[1], sex_ratio = sex_ratio, 
                                     pred_rate = pred_rate)
            

            # process_error(double z, double Y_m)
            error_list <- process_error(X_t1 = xtr, Y_t1 = ytr,
                                        X_t = X_0r[,i], Y_t = Y_0r[,i],
                                        sigma_x = sigma_x, sigma_y = sigma_y, rho = rho,
                                        z = z, Y_m = Y_0r[nrow(Y_0r),i],
                                        total_stages = n_aphid_stages + n_wasp_stages,
                                        living_aphids = n_aphid_stages + mum_days[1],
                                        sigma_d_mult = sigma_d_mult)
            xtr <- error_list$aphids
            ytr <- error_list$wasps

            error_list <- process_error(X_t1 = xts, Y_t1 = yts,
                                        X_t = X_0s[,i], Y_t = Y_0s[,i],
                                        sigma_x = sigma_x, sigma_y = sigma_y, rho = rho,
                                        z = z, Y_m = Y_0s[nrow(Y_0s),i],
                                        total_stages = n_aphid_stages + n_wasp_stages,
                                        living_aphids = n_aphid_stages + mum_days[1],
                                        sigma_d_mult = sigma_d_mult)
            xts <- error_list$aphids
            yts <- error_list$wasps
            
            # Do harvesting if it's in the list of harvest times
            # harvest()
            if (t %in% harvest_times[i,1:(ncol(harvest_times)-1)]) {
                # Kill non-parasitized aphids
                xtr <- harvest_surv * xtr
                xts <- harvest_surv * xts
                # Kill parasitized (but still living) aphids at the same rate
                ytr[1:mum_days[1]] <- harvest_surv * ytr[1:mum_days[1]]
                yts[1:mum_days[1]] <- harvest_surv * yts[1:mum_days[1]]
                # Kill all mummies
                ytr[(mum_days[1]+1):(length(ytr)-1)] <- 0
                yts[(mum_days[1]+1):(length(yts)-1)] <- 0
            }
            
            # Filling in values for the next iteration
            # (X_0r, X_0s, Y_0r, Y_0s hold time t info)
            X_0r[,i] <- xtr
            X_0s[,i] <- xts
            Y_0r[,i] <- ytr
            Y_0s[,i] <- yts
        }
        
        # Dispersal
        # ----
        
        # Aphids (first 4 instars apparently don't disperse)
        disp_stages <- (sum(instar_days[[clone[1,1]]][1:4])+1):nrow(X_0r)
        X_0r <- dispersal(X_0r, disp_aphid, disp_stages-1)  # -1 to make them C++ indices
        disp_stages <- (sum(instar_days[[clone[2,1]]][1:4])+1):nrow(X_0s)
        X_0s <- dispersal(X_0s, disp_aphid, disp_stages-1)
        
        # Wasps (only adults disperse)
        disp_stages <- nrow(Y_0r)
        Y_0r <- dispersal(Y_0r, disp_wasp, disp_stages-1)
        disp_stages <- nrow(Y_0s)
        Y_0s <- dispersal(Y_0s, disp_wasp, disp_stages-1)
        
        # X is sum of all living aphids
        Xr[t,] <- colSums(X_0r) + colSums(Y_0r[1:mum_days[1],])
        Xs[t,] <- colSums(X_0s) + colSums(Y_0s[1:mum_days[1],])
        
        # Y is proportion of live aphids that are parasitized
        Yr[t,] <- colSums(Y_0r[1:mum_days[1],]) / Xr[t,]
        Ys[t,] <- colSums(Y_0s[1:mum_days[1],]) / Xs[t,]
        
    }
    

    out_list <- list(Xr = Xr, Xs = Xs, Yr = Yr, Ys = Ys)

    return(out_list)
}




set.seed(60253704)
out_list <- HighTunnelExptSimfunct(X_0r,X_0s,Y_0r,Y_0s,a,attack_surv,K,K_y,k,h,s_y,
                                   # sigma_x,sigma_y,rho,1, # <- Full error
                                   # 0,0,0,1,  # <- No environmental error
                                   # sigma_x,sigma_y,rho,0, # <- No demographic error
                                   0,0,0,0, # <- No error
                                   max_time,n_fields,harvest_surv,harvest_times,disp_aphid,
                                   disp_wasp,pred_rate,
                                   n_aphid_stages, n_wasp_stages, instar_days, mum_days, 
                                   surv_juv, surv_adult, 
                                   rel_attack, sex_ratio, leslie_r, leslie_s, clone)
# Assigning out_list values to global ones for objects Xr, Xs, Yr, Ys
Xr <- out_list$Xr
Xs <- out_list$Xs
Yr <- out_list$Yr
Ys <- out_list$Ys

# base_p()

gg_p()





