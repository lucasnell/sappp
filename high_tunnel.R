# 
# This script seeks to replicate the high-tunnel model that Tony sent me Aug 2017
# 
# by:  Lucas A. Nell
# 







library(Matrix) # for sparse matrices


# Equivalent to diag(v, k) in MATLAB
diag_ml <- function(v, k = 0) {
    k <- as.integer(k)
    md <- diag(as.vector(v))
    md2 <- matrix(0, nrow(md) + abs(k), ncol(md) + abs(k))
    if (k > 0) {
        md2[1:nrow(md), (k+1):ncol(md2)] <- md
    } else if (k < 0) {
        md2[(abs(k)+1):nrow(md2), 1:ncol(md)] <- md
    } else {
        md2 <- md
    }
    return(md2)
}

# Check if any numbers are actually complex (imaginary part is != 0)
any_complex <- function(x) {
    any(sapply(x, function(xx) !identical(Im(xx), 0)))
}



# Expand from values per instar to per stage (i.e., day)
instar_to_stage <- function(stage_values, .n_stages = n_aphid_stages, 
                            .stage_days = stage_days) {
    
    .n_lines <- nrow(stage_days)
    
    .one_col <- function(jv, sdj, ns) c(rep(jv, sdj), rep(0, max(c(0, ns - sum(sdj)))))
    
    do.call(rbind, 
            lapply(1:.n_lines, function(j) {
                j_values <- stage_values[ifelse(nrow(stage_values) == 1, 1, j),]
                .one_col(j_values, .stage_days[j,], .n_stages)
            }))
    
}




# global n_aphid_stages stage_days mum_days surv_juv surv_adult rel_attack 
# sex_ratio Leslie_r Leslie_s clone


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




# =============================================
# lab parameters
# =============================================

n_aphid_stages <- 32
n_lines <- 2

n_wasp_stages <- sum(mum_days) + 1

# Number of days per instar and mummy development stage
stage_days <- rbind(c(2, 2, 2, 2, 19), c(1, 1, 1, 2, 23))
mum_days <- cbind(7, 3)


# juvenile survival
surv_juv <- cbind(0.9745, 0.9849)

# adult survival
surv_adult <- rbind(c(1.0000, 0.9949, 0.9818, 0.9534, 0.8805, 0.8367, 0.8532, 0.8786, 0.8823, 
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
rel_attack <- instar_to_stage(cbind(0.12, 0.27, 0.39, 0.16, 0.06))




sex_ratio <- 0.5


a <- 2.5        # parasitoid attack rate (2.32 in paper)
k <- 0.0005     # aphid density dependence (0.000467 in paper)
kp <- 0.0006    # parasitized aphid density dependence (0.00073319 in paper)
kk <- 0.1811    # ?? (difference between r at 20 and 27 deg??)
h <- 0.0363     # ?? (h was 0.008 in paper)

sw <- 0.55      # r at 27 degrees C ?? (0.554 in paper)

# used only in stochastic part:
# s1 <- 0
# s2 <- 0
# s3 <- 0
# rho <- 2 / (1 + exp(-s3)) - 1

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

# Create Leslie matrix from aphid info
Leslie_matrix <- function(n_stages, stage_days, clone_row, surv_juv, surv_adult, repro) {
    juv_time <- sum(stage_days[clone_row[1], 1:(ncol(stage_days)-1)])
    # Age-specific survivals
    LL <- diag_ml(c(surv_juv[,clone_row[1]] * matrix(1,1,juv_time), 
                    surv_adult[clone_row[2], 1:(n_stages-juv_time-1)]), -1)
    # Age-specific fecundities
    LL[1,(juv_time+1):(juv_time+stage_days[clone_row[1],ncol(stage_days)])] <- 
        repro[clone_row[2],1:(stage_days[clone_row[1],ncol(stage_days)])]
    return(LL)
}

# Not sure what this object creates yet, but it's used below
Leslie_dist <- function(L) {
    L_eigen <- eigen(L)
    SAD <- L_eigen$vectors
    r <- matrix(L_eigen$values, ncol = 1)
    rmax <- max(abs(r))

    SADdist <- SAD[, abs(r) == rmax]
    SADdist <- SADdist / sum(SADdist)
    if (!any_complex(SADdist)) SADdist <- as.numeric(SADdist)
    SADdist <- matrix(SADdist, ncol = 1)
    
    return(SADdist)
}

# resistant clones
# ------
Leslie_r <- Leslie_matrix(n_aphid_stages, stage_days, clone[1,], surv_juv, surv_adult, repro)
SADdistr <- Leslie_dist(Leslie_r)


# susceptible clones
# ------
Leslie_s <- Leslie_matrix(n_aphid_stages, stage_days, clone[2,], surv_juv, surv_adult, repro)
SADdists <- Leslie_dist(Leslie_s)

# Not sure what SADdist objects are


# =============================================
# field parameters: This is set up to have different harvesting patterns
# between n_fields fields.
# =============================================

# number of fields
n_fields <- 2

# kill rate at harvesting
kill <- 0.05

# dispersal rates between fields for aphids, adult wasps, and predators
disp_aphid <- 0.05
disp_wasp <- 1
disp_pred <- 0.8

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

init_x_resist <- prop_resist * init_x * SADdistr %*% matrix(1,1,n_fields)
init_x_susc <- (1-prop_resist) * init_x * SADdists %*% matrix(1,1,n_fields)

yr <- init_y * rbind(matrix(0, sum(mum_days), n_fields), c(1, 1))
ys <- init_y * rbind(matrix(0, sum(mum_days), n_fields), c(1, 1))



max_time <- cycle_length * (1 + n_cycles)
harvest_times <- rbind(c(cycle_length * 1:n_cycles),
                      c(cycle_length * 1, cycle_length * (2:(n_cycles-1)) - cycle_length/2,
                        cycle_length * n_cycles))

xs <- init_x_susc
xr <- init_x_resist










HighTunnelExptSimfunct <- function(
    xr,xs,yr,ys,a,resist_surv,k,kp,kk,h,sw,
    # s1,s2,rho,
    # sm,
    max_time,n_fields,kill,
    harvest_times,disp_aphid,disp_wasp,disp_pred,
    n_aphid_stages, n_wasp_stages, stage_days, mum_days, surv_juv, surv_adult, 
    rel_attack, sex_ratio, Leslie_r, Leslie_s, clone) {
    
    # global n_aphid_stages stage_days mum_days surv_juv surv_adult 
    # rel_attack sex_ratio Leslie_r Leslie_s clone

    # This is for measurement error part:
    # Su <- 0.2015^2 * diag(1, 2)
    # # increased ME for mummies for development time
    # Su[2,2] <- sm * Su[2,2]
    
    Xr <- matrix(0, max_time, n_fields)
    Yr <- matrix(0, max_time, n_fields)
    Xs <- matrix(0, max_time, n_fields)
    Ys <- matrix(0, max_time, n_fields)
    
    LLr <- Matrix(Leslie_r, sparse = TRUE)
    LLs <- Matrix(Leslie_s, sparse = TRUE)
    
    for (t in 1:max_time) {
        for (i in 1:n_fields) {
            ypr <- yr[1:mum_days[1], i]
            ypr <- ypr[ypr>0]
            yps <- ys[1:mum_days[1], i]
            yps <- yps[yps>0]
            
            sxyt <- sum(c(xr[,i], xs[,i], ypr, yps))
            Kt <- 1 / (1 + k * sxyt)
            Kpt <- 1 / (1 + kp * sxyt)
            
            
            mm <- rbind(a * rel_attack[clone[1,1],] * (yr[nrow(yr),i] + ys[nrow(ys),i]) / 
                            (h * sxyt + 1))
            AA <- (1 + t(mm) / kk)
            As <- AA^(-kk)
            
            mm <- rbind(a * rel_attack[clone[2,1],] * (yr[nrow(yr),i] + ys[nrow(ys),i]) / 
                            (h * sxyt + 1))
            AA <- (1 + t(mm) / kk)
            Ar <- AA^(-kk) + resist_surv[1] * t(mm) * AA^(-kk-1) + resist_surv[2] *
                (1-(AA^(-kk) + t(mm) * AA^(-kk-1)))
            
            xtr <- (disp_pred * Kt * Ar) * as.matrix(LLr %*% xr[,i])
            xts <- (disp_pred * Kt * As) * as.matrix(LLs %*% xs[,i])
            
            yt <- cbind(yr[,i])
            y <- cbind(yr[,i])
            x <- cbind(xr[,i])
            yt[length(yt)] <- sw * y[length(y)] + sex_ratio * y[(length(y)-1)]
            yt[(mum_days[1]+2):(length(yt)-1)] <- disp_pred * y[(mum_days[1]+1):(length(y)-2)]
            yt[2:(mum_days[1]+1)] <- Kpt * surv_juv[clone[1,1]] * y[1:mum_days[1]]
            yt[1] <- (Kpt * t(matrix(1,n_aphid_stages,1) - Ar)) %*% as.matrix(LLr %*% x)
            ytr <- yt
            
            yt <- cbind(ys[,i])
            y <- cbind(ys[,i])
            x <- cbind(xs[,i])
            yt[length(yt)] <- sw * y[length(y)] + sex_ratio * y[(length(y)-1)]
            yt[(mum_days[1]+2):(length(yt)-1)] <- disp_pred * y[(mum_days[1]+1):(length(y)-2)]
            yt[2:(mum_days[1]+1)] <- Kpt * surv_juv[clone[2,1]] * y[1:mum_days[1]]
            yt[1] <- (Kpt * t(matrix(1,n_aphid_stages,1) - As)) %*% as.matrix(LLs %*% x)
            yts <- yt
            
            # This code for stochasticity is dead
            # process error for aphids, parasitized aphids, and adult parasitoids		
            # Se <- matrix(0, n_aphid_stages+n_wasp_stages, n_aphid_stages+n_wasp_stages)
            # Se[1:nap,1:nap] <- 
            # 	s1^2 * (rho * matrix(1,nap,nap) + (1-rho) * diag(1,nap))
            # 
            # Se[(nap+1):(end-1),(nap+1):(end-1)] <- 0
            # 
            # Se[1:nap,(nap+1):(end-1)] <- 0
            # Se[(nap+1):(end-1),1:nap] <- 0
            # 
            # Se[nrow(Se),ncol(Se)] <- s3^2
            # ypick <- which(diag(Se) > 0)
            # 
            # iD <- t(chol(Se[ypick,ypick]))
            # E <- iD %*% rnorm(length(ypick))
            
            if (t %in% harvest_times[i,1:(ncol(harvest_times)-1)]) {
                # Tony used a matrix for the end of a sequence (1:mum_days).
                # The behavior of this in matlab is to just end at the first element of 
                # the matrix (i.e., 1:mum_days[1,1])
                mumInts_ <- mum_days[1,1]
                xtr <- kill * xtr
                ytr[1:mumInts_] <- kill * ytr[1:mumInts_]
                ytr[(mumInts_+1):(length(ytr)-1)] <- 0
                
                xts <- kill * xts
                yts[1:mumInts_] <- kill * yts[1:mumInts_]
                yts[(mumInts_+1):(length(yts)-1)] <- 0
            }
            
            xr[,i] <- xtr
            xs[,i] <- xts
            yr[,i] <- ytr
            ys[,i] <- yts
        }
        
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
        
        nap <- n_aphid_stages + mum_days[1]
        
        yy <- rbind(xr, yr)
        Xr[t,] <- colSums(yy[1:nap,])
        Yr[t,] <- colSums(yy[(n_aphid_stages+1):nap,]) / colSums(yy[1:nap,])
        
        yy <- rbind(xs, ys)
        Xs[t,] <- colSums(yy[1:nap,])
        Ys[t,] <- colSums(yy[(n_aphid_stages+1):nap,]) / colSums(yy[1:nap,])
    }
    

    out_list <- list(Xr = Xr, Xs = Xs, Yr = Yr, Ys = Ys, xr = xr, xs = xs, 
                     yr = yr, ys = ys)
    
    return(out_list)
}



out_list <- HighTunnelExptSimfunct(xr,xs,yr,ys,a,resist_surv,k,kp,kk,h,sw,
                                   # s1,s2,rho,
                                   # sm,
                                   max_time,n_fields,kill,harvest_times,disp_aphid,
                                   disp_wasp,disp_pred,
                                   n_aphid_stages, n_wasp_stages, stage_days, mum_days, 
                                   surv_juv, surv_adult, 
                                   rel_attack, sex_ratio, Leslie_r, Leslie_s, clone)
# Assigning out_list values to global ones for objects Xr, Xs, Yr, Ys, xr, xs, yr, & ys
invisible(
    lapply(names(out_list), 
           function(NAME) {
               eval(parse(text = paste0(NAME, ' <<- out_list$', NAME)))
           }))


Ymax <- 1.2 * max(Xr+Xs)
Ymax2 <- 1.2 * max(Xr+Xs)
range <- c(0, max_time, 0, Ymax)
min_time <- 1

# Figure 1
plot(1:(max_time-min_time+1),Xr[min_time:max_time,1]+Xs[min_time:max_time,1], 
     type = 'l', col = 'dodgerblue')
lines(1:(max_time-min_time+1),Xr[min_time:max_time,2]+Xs[min_time:max_time,2], 
      col = 'dodgerblue', lty = 2)
lines(1:(max_time-min_time+1),Ymax2 * (Yr[min_time:max_time,1]+Ys[min_time:max_time,1]), 
      col = 'firebrick')
lines(1:(max_time-min_time+1),Ymax2 * (Yr[min_time:max_time,2]+Ys[min_time:max_time,2]), 
      col = 'firebrick', lty = 2)
lines(1:(max_time-min_time+1),Ymax2*(Xr[min_time:max_time,1]+Xr[min_time:max_time,2])/
         (Xr[min_time:max_time,1]+Xr[min_time:max_time,2]+Xs[min_time:max_time,1]+
              Xs[min_time:max_time,2]), 
      col = 'black')

# Blue is aphid abundances
# Red is parasitoid abundances
# Black is (resistant aphids) / (total aphids)
# Dotted lines might be a different harvesting regime

