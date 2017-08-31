# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# # library(Matrix) # for sparse matrices
# 
# # Provides functions instar_to_stage, leslie_matrix, leslie_sad, attack_probs, 
# # parasitoid_abunds
Rcpp::sourceCpp('high_tunnel.cpp')

source('classes.R')




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

clone <- data.frame(juv = c(2, 2), adult = c(1, 2))
rownames(clone) <- c('resistant', 'susceptible')

# Values dictate which rows will be selected from matrices of demographic rates below
# (instar_days, surv_juv, surv_adult, repro)



# =============================================
# lab parameters
# =============================================

n_aphid_stages <- 32
n_lines <- 2


# Number of days per instar (for low and high growth rates)
instar_days <- list(low = as.integer(c(2, 2, 2, 2, 19)), 
                    high = as.integer(c(1, 1, 1, 2, 23)))
# Number of days for parasitized (but still living) aphids and mummies
mum_days <- cbind(7, 3)

n_wasp_stages <- sum(mum_days) + 1

# juvenile survival
surv_juv <- list(low = 0.9745, high = 0.9849)

# adult survival
surv_adult <- list(
    low = rbind(c(1.0000, 0.9949, 0.9818, 0.9534, 0.8805, 0.8367, 0.8532, 0.8786, 
                  0.8823, 0.8748, 0.8636, 0.8394, 0.8118, 0.8096, 0.8240, 0.8333, 
                  0.7544, 0.5859, 0.4155, 0.2216, rep(0, 180))),
    high = rbind(c(1.0000, 0.9986, 0.9951, 0.9874, 0.9675, 0.9552, 0.9550, 0.9549, 
                   0.9462, 0.8992, 0.8571, 0.8408, 0.8281, 0.8062, 0.7699, 0.7500, 
                   0.7559, 0.7649, 0.7240, 0.4367, rep(0, 180))))

# reproduction
repro <- list(
    low = rbind(c(0, 2.5925, 4.4312, 5.1403, 5.5190, 5.6633, 5.6010, 5.4577, 5.2904, 
                  5.0613, 4.6970, 3.3577, 1.5946, 1.0817, 0.9666, 0.8333, 0.4689, 
                  0.0709, 0, 0, 0, 0, rep(0, 178))),
    high = rbind(c(0, 3.1975, 5.4563, 6.2996, 6.7372, 6.9030, 6.8210, 6.6100, 
                   6.1962, 5.1653, 4.1837, 3.6029, 3.1023, 2.4799, 1.6909, 1.1750, 
                   1.0148, 0.9096, 0.7821, 0.6430, 0.5000, 0.3531, rep(0, 178))))


# Relative attack rate on the different instars from Ives et al 1999
# `instar_to_stage` converts these values from per-instar to per-day
rel_attack <- list(low = instar_to_stage(rbind(0.12, 0.27, 0.39, 0.16, 0.06), 
                                         n_aphid_stages, instar_days$low),
                   high = instar_to_stage(rbind(0.12, 0.27, 0.39, 0.16, 0.06),
                                          n_aphid_stages, instar_days$high))


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


# These are the survivals of singly attacked and multiply attacked
# resistant aphids
attack_surv <- cbind(0.9, 0.6)



# =============================================
# set up Leslie matrices
# =============================================

# resistant clones
# ------
leslie_r <- leslie_matrix(n_aphid_stages, instar_days[[clone[1,1]]],
                          surv_juv[[clone[1,1]]],
                          surv_adult[[clone[1,2]]], repro[[clone[1,2]]])
sad_r <- leslie_sad(leslie_r)


# susceptible clones
# ------
leslie_s <- leslie_matrix(n_aphid_stages, instar_days[[clone[2,1]]],
                          surv_juv[[clone[2,1]]],
                          surv_adult[[clone[2,2]]], repro[[clone[2,2]]])
sad_s <- leslie_sad(leslie_s)



# =============================================
# field parameters: This is set up to have different harvesting patterns
# between n_fields fields.
# =============================================

# number of fields
n_fields <- 2

# survival rate at harvesting
harvest_surv <- 0.05

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
X_0r <- prop_resist * init_x * sad_r %*% matrix(1,1,n_fields)
X_0s <- (1-prop_resist) * init_x * sad_s %*% matrix(1,1,n_fields)

# Initial parasitoid densities by stage (starting with no parasitized aphids or mummies)
Y_0r <- init_y * rbind(matrix(0, sum(mum_days), n_fields), c(1, 1))
Y_0s <- init_y * rbind(matrix(0, sum(mum_days), n_fields), c(1, 1))

# Setting total time (days) and times for harvesting
max_time <- cycle_length * (1 + n_cycles)
harvest_times <- rbind(c(cycle_length * 1:n_cycles),
                       c(cycle_length * 1, cycle_length * (2:(n_cycles-1)) - cycle_length/2,
                         cycle_length * n_cycles))



# leslie_r <- leslie_matrix(n_aphid_stages, instar_days$high,
#                           surv_juv$high, surv_adult$high, repro$low)
# prop_resist * init_x * leslie_sad(leslie_r)


# =====================================================================================
# =====================================================================================

#       APHID LINES

# =====================================================================================
# =====================================================================================

# -------
# Resistant aphid line info
# -------
res_line <- aphid_const$new(
    leslie = leslie_matrix(instar_days$high, surv_juv$high, surv_adult$high, repro$low),
    rel_attack = rel_attack$high,
    a = a, K = K, K_y = K_y, k = k, h = h, s_y = s_y, sigma_x = sigma_x, 
    sigma_y = sigma_y, rho = rho,
    disp_stages = (sum(instar_days[[clone[2,1]]][1:4])+1):32, 
    mum_days = mum_days, 
    aphid_density_0 = prop_resist * init_x,
    attack_surv = attack_surv)
res_line




# -------
# Susceptible aphid line info
# -------
# Difference is no resistance, higher reproduction, higher starting density
sus_line <- aphid_const$new(
    leslie = leslie_matrix(instar_days$high, surv_juv$high, surv_adult$high, repro$high),
    rel_attack = rel_attack$high,
    a = a, K = K, K_y = K_y, k = k, h = h, s_y = s_y, sigma_x = sigma_x, 
    sigma_y = sigma_y, rho = rho,
    disp_stages = (sum(instar_days[[clone[2,1]]][1:4])+1):32, 
    mum_days = mum_days, 
    aphid_density_0 = (1 - prop_resist) * init_x)
sus_line






# =====================================================================================
# =====================================================================================

#       FIELDS

# =====================================================================================
# =====================================================================================







# ====================================================================================
# ====================================================================================

# Plotting

# ====================================================================================
# ====================================================================================


base_p <- function(ymult = 1) {
    Ymax <- ymult * max(Xr[,1]+Xs[,1])
    min_time <- 1
    
    # Figure 1
    plot(1:(max_time-min_time+1),Xr[min_time:max_time,1]+Xs[min_time:max_time,1],
         type = 'l', col = 'dodgerblue', ylab = '', xlab = 'time', main = '')
    lines(1:(max_time-min_time+1),Xr[min_time:max_time,2]+Xs[min_time:max_time,2],
          col = 'dodgerblue', lty = 2)
    lines(1:(max_time-min_time+1),Ymax * (Yr[min_time:max_time,1]+Ys[min_time:max_time,1]),
          col = 'firebrick')
    lines(1:(max_time-min_time+1),Ymax * (Yr[min_time:max_time,2]+Ys[min_time:max_time,2]),
          col = 'firebrick', lty = 2)
    lines(1:(max_time-min_time+1),Ymax*rowSums(Xr[min_time:max_time,])/
              (rowSums(Xr[min_time:max_time,]) + rowSums(Xs[min_time:max_time,])),
          col = 'black')
    
    # Blue is aphid abundances
    # Red is parasitoid abundances
    # Black is (resistant aphids) / (total aphids)
    # Dotted lines are field # 2 (different harvesting regime)
}




gg_p <- function() {
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
}