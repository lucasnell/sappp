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

# source('parameters.R')




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











HighTunnelExptSimfunct <- function(
    xr,xs,yr,ys,a,resist_surv,K,K_y,k,h,s_y,
    sigma_x,sigma_y,rho,sigma_d_mult,
    # sm,
    max_time,n_fields,kill,
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
    LLr <- leslie_r  # Matrix(leslie_r, sparse = TRUE)
    LLs <- leslie_s  # Matrix(leslie_s, sparse = TRUE)
    
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
            
            # X(t+1) (first line of equation 2; pred_rate was added after paper)
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
            
            # Do harvesting if it's in the list of harvest times
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
        
        # Dispersal
        # ----
        
        # Aphids (first 4 instars apparently don't disperse)
        disp_stages <- (sum(instar_days[clone[1,1],1:4])+1):nrow(xr)
        xr <- dispersal(xr, disp_aphid, disp_stages-1)  # -1 to make them C++ indices
        disp_stages <- (sum(instar_days[clone[2,1],1:4])+1):nrow(xs)
        xs <- dispersal(xs, disp_aphid, disp_stages-1)
        
        # Wasps (only adults disperse)
        disp_stages <- nrow(yr)
        yr <- dispersal(yr, disp_wasp, disp_stages-1)
        disp_stages <- nrow(ys)
        ys <- dispersal(ys, disp_wasp, disp_stages-1)
        
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
                                   n_aphid_stages, n_wasp_stages, instar_days, mum_days, 
                                   surv_juv, surv_adult, 
                                   rel_attack, sex_ratio, leslie_r, leslie_s, clone)
# Assigning out_list values to global ones for objects Xr, Xs, Yr, Ys
Xr <- out_list$Xr
Xs <- out_list$Xs
Yr <- out_list$Yr
Ys <- out_list$Ys

# base_plot()




aphids <- as_data_frame(cbind(Xr, Xs)) %>%
    mutate(time = 1:nrow(Xs), `1` = V1+V3, `2` = V2+V4,
           r_prop = (V1+V2)/(V1+V2+V3+V4)) %>% # <-- (resistant aphids) / (total aphids)
    select(-1:-4) %>% 
    gather('field', 'density', 2:3, factor_key = TRUE)

wasps <- as_data_frame(cbind(Ys, Yr)) %>%
    mutate(time = 1:nrow(Ys), `1` = (V1+V3)/2, `2` = (V2+V4)/2) %>% 
    select(-1:-4) %>% 
    gather('field', 'density', 2:3, factor_key = TRUE)

axis_mult = max(aphids$density)
# axis_mult = 330

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
