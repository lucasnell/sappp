
# Constant info about an aphid line (and wasps that parasitize them)
# This should be used among fields for the same aphid line
aphid_const <- setRefClass(
    
    'aphid_const', 
    
    fields = list(
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
        mum_days = 'integer'        # number of days per mummy stage (aphid alive & dead)
    )
)

# Info about an aphid line (and wasps that parasitize them), some of which changes 
# through time.
# This should be used for a single line in a single field.
aphid_iter <- setRefClass(
    
    'aphid_iter',
    
    fields = list(
        line_info = 'aphid_const',  # Aphid line info (doesn't change through time)
        X_t = 'numeric',            # Aphid density at time t
        X_t1 = 'numeric',           # Aphid density at time t+1
        Y_t = 'numeric',            # Wasp density at time t
        Y_t1 = 'numeric',           # Wasp density at time t+1
        A = 'numeric'               # Attack probabilities at time t+1
    )
)


# Info about a patch that's constant through time
patch_const <- setRefClass(
    
    'patch_const', 
    
    fields = list(
        harvest_times = 'integer',  # times when harvesting occurs
        max_time = 'integer'        # number of time points to simulate
    )
)

# Info about a patch (and its inhabitants) that changes through time
patch_iter <- setRefClass(
    
    'patch_iter',
    
    field = list(
        patch_info = 'patch_const', # Patch info (doesn't change through time)
        z = 'numeric',              # Sum of all living aphids at time t
        x = 'numeric',              # Sum of non-parasitized aphids at time t
        S = 'numeric',              # Density dependence for aphids
        S_y = 'numeric'             # Density dependence for wasps
    )
    
)


