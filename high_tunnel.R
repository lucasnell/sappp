# 
# This script seeks to replicate the high-tunnel model that Tony sent me Aug 2017
# 
# by:  Lucas A. Nell
# 


library(Matrix) # for sparse matrices

# Provides functions instar_to_stage, leslie_matrix, leslie_sad, attack_probs, 
# parasitoid_abunds
Rcpp::sourceCpp('high_tunnel.cpp')




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
k <- 0.0005     # aphid density dependence (K; 0.000467)
kp <- 0.0006    # parasitized aphid density dependence (K_y; 0.00073319)
kk <- 0.1811    # aggregation parameter of the negative binomial distribution (k; 0.35)
h <- 0.0363     # parasitoid attack rate handling time (h; 0.008, 0.029 at 27 deg C)
sw <- 0.55      # parasitoid adult daily survival (s_y; 0.69)
s1 <- 0.44  # 0  # environmental standard deviation for aphids (sigma_x; 0.44)
# s2 <- 0  # Doesn't appear to do anything
s3 <- 0.35  # 0  # environmental standard deviation for parasitoids (sigma_y; 0.70)
rho <- 2 / (1 + exp(-s3)) - 1  # environmental correlation among instars (rho; 1.0)




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
    xr,xs,yr,ys,a,resist_surv,k,kp,kk,h,sw,
    # s1,s2,rho,
    # sm,
    max_time,n_fields,kill,
    harvest_times,disp_aphid,disp_wasp,pred_rate,
    n_aphid_stages, n_wasp_stages, stage_days, mum_days, surv_juv, surv_adult, 
    rel_attack, sex_ratio, leslie_r, leslie_s, clone) {
    
    # This is for measurement error part:
    # Su <- 0.2015^2 * diag(1, 2)
    # # increased ME for mummies for development time
    # Su[2,2] <- sm * Su[2,2]
    
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
            # Numbers of parasitized (but still living) aphids
            ypr <- yr[1:mum_days[1], i]
            ypr <- ypr[ypr>0]
            yps <- ys[1:mum_days[1], i]
            yps <- yps[yps>0]
            
            # Total number of living aphids (z in the paper)
            z <- sum(c(xr[,i], xs[,i], ypr, yps))
            
            # Equivalent to S(z(t)) and S_y(z(t)) [no idea why equations are different]
            Kt <- 1 / (1 + k * z)
            Kpt <- 1 / (1 + kp * z)
            
            # Matrices of attack probabilities using equation 6 from paper
            As <- attack_probs(a = a, p_i = rel_attack[clone[1,1],], 
                               Y_m = yr[nrow(yr),i] + ys[nrow(ys),i], 
                               x = z, h = h, k = kk, resist_surv = numeric(0))
            Ar <- attack_probs(a = a, p_i = rel_attack[clone[2,1],], 
                               Y_m = yr[nrow(yr),i] + ys[nrow(ys),i], 
                               x = z, h = h, k = kk, resist_surv = resist_surv)
            
            # X(t+1) (first line of equation 2; pred_rate was added later)
            xtr <- (pred_rate * Kt * Ar) * as.matrix(LLr %*% xr[,i])
            xts <- (pred_rate * Kt * As) * as.matrix(LLs %*% xs[,i])
            
            # Filling in column of parasitoid stage abundances
            ytr <- parasitoid_abunds(S_y_zt = Kpt, A = Ar, L = LLr, X = xr[,i], 
                                     Y_t = yr[,i], s_i = surv_juv[clone[1,1]], 
                                     s_y = sw, m_1 = mum_days[1], sex_ratio = sex_ratio, 
                                     pred_rate = pred_rate)
            yts <- parasitoid_abunds(S_y_zt = Kpt, A = As, L = LLs, X = xs[,i], 
                                     Y_t = ys[,i], s_i = surv_juv[clone[2,1]], 
                                     s_y = sw, m_1 = mum_days[1], sex_ratio = sex_ratio, 
                                     pred_rate = pred_rate)
            
            
            # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # # This code for stochasticity is dead
            # # process error for aphids, parasitized aphids, and adult parasitoids
            # nap <- n_aphid_stages + mum_days[1]
            # Se <- matrix(0, n_aphid_stages+n_wasp_stages, n_aphid_stages+n_wasp_stages)
            # Se[1:nap,1:nap] <-
            # 	# s1^2 * (rho * matrix(1,nap,nap) + (1-rho) * diag(1,nap))
            #     (s1^2 + min(0.5, 1 / abs(1 + sum(yy(1:nap))))) * 
            #     (rho * matrix(1,nap,nap)+(1-rho)*diag(1,nap));
            # 
            # Se[(nap+1):(nrow(Se)-1),(nap+1):(ncol(Se)-1)] <- 0
            # 
            # Se[1:nap,(nap+1):(ncol(Se)-1)] <- 0
            # Se[(nap+1):(nrow(Se)-1),1:nap] <- 0
            # 
            # Se[nrow(Se),ncol(Se)] <- s3^2
            # ypick <- which(diag(Se) > 0)
            # 
            # iD <- t(chol(Se[ypick,ypick]))
            # E <- iD %*% rnorm(length(ypick))
            # 
            # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
            
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
    

    out_list <- list(Xr = Xr, Xs = Xs, Yr = Yr, Ys = Ys, xr = xr, xs = xs, 
                     yr = yr, ys = ys)
    
    return(out_list)
}



out_list <- HighTunnelExptSimfunct(xr,xs,yr,ys,a,resist_surv,k,kp,kk,h,sw,
                                   # s1,s2,rho,
                                   # sm,
                                   max_time,n_fields,kill,harvest_times,disp_aphid,
                                   disp_wasp,pred_rate,
                                   n_aphid_stages, n_wasp_stages, stage_days, mum_days, 
                                   surv_juv, surv_adult, 
                                   rel_attack, sex_ratio, leslie_r, leslie_s, clone)
# Assigning out_list values to global ones for objects Xr, Xs, Yr, Ys, xr, xs, yr, & ys
invisible(
    lapply(names(out_list), 
           function(NAME) {
               eval(parse(text = paste0(NAME, ' <<- out_list$', NAME)))
           }))




library(dplyr)
library(tidyr)
library(ggplot2)


aphids <- as_data_frame(cbind(Xr, Xs)) %>%
    mutate(time = 1:nrow(Xs), `1` = V1+V3, `2` = V2+V4,
           r_prop = (V1+V2)/(V1+V2+V3+V4)) %>% # <-- (resistant aphids) / (total aphids)
    select(-1:-4) %>% 
    gather('field', 'density', 2:3, factor_key = TRUE)

wasps <- as_data_frame(cbind(Ys, Yr)) %>%
    mutate(time = 1:nrow(Ys), `1` = V1+V3, `2` = V2+V4) %>% 
    select(-1:-4) %>% 
    gather('field', 'density', 2:3, factor_key = TRUE)


ggplot(aphids, aes(time)) +
    theme_classic() +
    theme(legend.position = c(0.01, 0.99), legend.direction = 'horizontal',
          legend.justification = c('left', 'top')) +
    geom_line(aes(y = density, linetype = field), color = 'dodgerblue') +
    geom_line(data = wasps, aes(y = density * 330, linetype = field), 
              color = 'firebrick') +
    geom_line(aes(y = r_prop * 330)) +
    geom_text(data = data_frame(time = c(350, 50, 300), y = c(185, 165, 230), 
                                lab = c('aphids', '% parasit.', '% resist.')),
              aes(y = y, label = lab), color = c('dodgerblue', 'firebrick', 'black'),
              hjust = 0, vjust = 0) +
    scale_y_continuous('aphid density', 
                       sec.axis = sec_axis(~ . * 100 / 330, name = '%'))






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

