

# Info about an aphid line (and wasps that parasitize them)
aphid_line <- setRefClass(
    
    'aphid_line', 
    
    fields = list(
        # --------
        # Constant info
        # --------
        leslie = 'matrix',          # Leslie matrix with survival and reproduction
        rel_attack = 'numeric',     # relative wasp attack rates by aphid stage
        X_0 = 'numeric',            # initial aphid abundances by stage
        Y_0 = 'numeric',            # initial wasp abundances by stage
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
        disp_wasp = 'numeric',      # dispersal rate for wasps
        pred_rate = 'numeric',      # predation on aphids and mummies
        n_aphid_stages = 'integer', # number of aphid stages (i.e., days)
        n_wasp_stages = 'integer',  # number of wasp stages (i.e., days)
        instar_days = 'integer',    # number of days per aphid instar
        mum_days = 'integer',       # number of days per mummy stage (aphid alive & dead)
        
        # --------
        # Changing info
        # --------
        X_t = 'numeric',            # Aphid density at time t
        X_t1 = 'numeric',           # Aphid density at time t+1
        Y_t = 'numeric',            # Wasp density at time t
        Y_t1 = 'numeric',           # Wasp density at time t+1
        A = 'numeric'               # Attack probabilities at time t+1
    )
)

# Info about a patch
patch <- setRefClass(
    
    'patch', 
    
    fields = list(
        # --------
        # Constant info
        # --------
        harvest_times = 'integer',  # times when harvesting occurs
        max_time = 'integer',       # number of time points to simulate
        
        # --------
        # Changing info
        # --------
        z = 'numeric',              # Sum of all living aphids at time t
        S = 'numeric',              # Density dependence for aphids
        S_y = 'numeric'             # Density dependence for wasps
    )
)








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
# (instar_days, surv_juv, surv_adult, repro)



# =============================================
# lab parameters
# =============================================

n_aphid_stages <- 32
n_lines <- 2


# Number of days per instar
instar_days <- rbind(c(2, 2, 2, 2, 19), c(1, 1, 1, 2, 23))
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
                              n_aphid_stages, instar_days)




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
resist_surv <- cbind(0.9, 0.6)




# =============================================
# set up Leslie matrices
# =============================================

# resistant clones
# ------
leslie_r <- leslie_matrix(n_aphid_stages, instar_days[clone[1,1],], surv_juv[clone[1,1]],
                          surv_adult[clone[1,2],], repro[clone[1,2],])
sad_r <- leslie_sad(leslie_r)


# susceptible clones
# ------
leslie_s <- leslie_matrix(n_aphid_stages, instar_days[clone[2,1],], surv_juv[clone[2,1]],
                          surv_adult[clone[2,2],], repro[clone[2,2],])

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

