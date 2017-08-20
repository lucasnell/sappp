# 
# This script seeks to replicate the high-tunnel model that Tony sent me Aug 2017
# 
# by:  Lucas A. Nell
# 


library(Matrix) # for sparse matrices

HighTunnelExptSimfunct_17Feb16 <- function(
    xr,xs,yr,ys,a,fresist,k,kp,kk,h,sw,s1,s2,rho,sm,Tmax,nfield,kill,
    harvesttimes,da,dw,pred) {
    
    # global totStage stageInts mumInts sj sa relatt sexratio Lr Ls clone

    Su <- 0.2015^2 * diag(1, 2)
    
    # increased ME for mummies for development time
    Su[2,2] <- sm * Su[2,2]
    
    nx  <- totStage
    ny <- sum(mumInts) + 1
    
    Xr <- matrix(0, Tmax, nfield)
    Yr <- matrix(0, Tmax, nfield)
    Xs <- matrix(0, Tmax, nfield)
    Ys <- matrix(0, Tmax, nfield)
    
    LLr <- Matrix(Lr, sparse = TRUE)
    LLs <- Matrix(Ls, sparse = TRUE)
    
    for (t in 1:Tmax) {
        for (i in 1:nfield) {
            ypr <- yr[1:mumInts[1], i]
            ypr <- ypr[ypr>0]
            yps <- ys[1:mumInts[1], i]
            yps <- yps[yps>0]
            
            sxyt <- sum(c(xr[,i], xs[,i], ypr, yps))
            Kt <- 1 / (1 + k * sxyt)
            Kpt <- 1 / (1 + kp * sxyt)
            
            
            mm <- a * relatt[clone[1,1],] * (yr[end,i] + ys[end, i]) / (h * sxyt + 1)
            AA <- (1 + t(mm) / kk)
            As <- AA^(-kk)
            
            mm <- a * relatt[clone[2,1],] * (yr[end,i] + ys[end,i]) / (h * sxyt + 1)
            AA <- (1 + t(mm) / kk)
            Ar <- AA^(-kk) + fresist[1] * t(mm) * AA^(-kk-1) + fresist[2] * 
                (1-(AA^(-kk) + t(mm) * AA^(-kk-1)))
            
            xtr <- pred * Kt * Ar * (LLr * xr[,i])
            xts <- pred * Kt * As * (LLs * xs[,i])
            
            yt <- yr[,i]
            y <- yr[,i]
            x <- xr[,i]
            yt[end] <- sw * y[end] + sexratio * y[end-1]
            yt[(mumInts[1]+2):(end-1)] <- pred * y[(mumInts[1]+1):(end-2)]
            yt[2:(mumInts[1]+1)] <- Kpt * sj[clone[1,1]] * y[1:mumInts[1]]
            yt[1] <- Kpt * t(matrix(1,totStage,1) - Ar) * (LLr * x)
            ytr <- yt
            
            yt <- ys[,i]
            y <- ys[,i]
            x <- xs[,i]
            yt[end] <- sw * y[end] + sexratio * y[end-1]
            yt[(mumInts[1]+2):(end-1)] <- pred * y[(mumInts[1]+1):(end-2)]
            yt[2:(mumInts[1]+1)] <- Kpt * sj[clone[2,1]] * y[1:mumInts[1]]
            yt[1] <- Kpt * t(matrix(1,totStage,1) - As) * (LLs * x)
            yts <- yt
            
            if (t %in% harvesttimes[i,1:(end-1)]) {
                xtr <- kill * xtr
                ytr[1:mumInts] <- kill * ytr[1:mumInts]
                ytr[(mumInts+1):(end-1)] <- 0
                
                xts <- kill * xts
                yts[1:mumInts] <- kill * yts[1:mumInts]
                yts[(mumInts+1):(end-1)] <- 0
            }
            
            xr[,i] <- xtr
            xs[,i] <- xts
            yr[,i] <- ytr
            ys[,i] <- yts
        }
        
        
    }
    
    dispersing <- da * t(mean(t(xr[(sum(stageInts[clone[1,1],1:4])+1):end,])))
    xr[(sum(stageInts[clone[1,1],1:4])+1):end,] <- (1-da) * 
        xr[(sum(stageInts[clone[1,1],1:4])+1):end,] + 
        dispersing * matrix(1,1,nfield)
    
    dispersing <- da * t(mean(t(xr[(sum(stageInts[clone[2,1],1:4])+1):end,])))
    xs[(sum(stageInts[clone[2,1],1:4])+1):end,] <- (1-da) * 
        xs[(sum(stageInts[clone[2,1],1:4])+1):end,] + 
        dispersing * matrix(1,1,nfield)
    
    dispersingw <- dw * mean(yr[end,])
    yr[end,] <- ((1-dw) * yr[end,] + dispersingw)
    
    dispersingw <- dw * mean(ys[end,])
    ys[end,] <- ((1-dw) * ys[end,] + dispersingw)
    
    nap <- nx + mumInts[1]
    
    yy <- c(xr, yr)
    Xr[t,] <- sum(yy[1:nap,]);
    Yr[t,] <- sum(yy[(nx+1):nap,]) / sum(yy[1:nap,])
    
    yy <- c(xs, ys)
    Xs[t,] <- sum(yy[1:nap,])
    Ys[t,] <- sum(yy[(nx+1):nap,]) / sum(yy[1:nap,])

}












# global totStage stageInts mumInts sj sa relatt sexratio Lr Ls clone


# % NOTE: The code is set up for 2 clones that have different life histories.
# % The life histories are taken from lab experiments and correspond to
# % demographic parameters at 20 and 27 C. See Meisner et al. 2014.
# 
# % The juvenile and adult demographic rates are separated, so you can pick
# % them separately for the different clones. "Resistant" clones are
# % resistant to parasitism.
# 
# % rows: 1=resistant, susceptible = 2
# % columns: 1=aphid juv, 2=aphid adult
# % values 1=low growth rate, 2=high growth rate
clone <- matrix(c(2, 1, 2, 2), 2, 2, byrow = TRUE)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % lab parameters

totStage <- 32

# % this sets the time scale in terms of the numbers of days per instar and
# % mummy development time
stageInts <- matrix(c(2, 2, 2, 2, 19,
                      1, 1, 1, 2, 23), 2, 5, byrow = TRUE)
mumInts <- c(7, 3)

# % juvenile survival
sj <- c(0.9745, 0.9849)

# % adult survival
sa=rbind(c(1.0000, 0.9949, 0.9818, 0.9534, 0.8805, 0.8367, 0.8532, 0.8786, 0.8823, 
           0.8748, 0.8636, 0.8394, 0.8118, 0.8096, 0.8240, 0.8333, 0.7544, 0.5859, 
           0.4155, 0.2216),
         c(1.0000, 0.9986, 0.9951, 0.9874, 0.9675, 0.9552, 0.9550, 0.9549, 0.9462, 
           0.8992, 0.8571, 0.8408, 0.8281, 0.8062, 0.7699, 0.7500, 0.7559, 0.7649, 
           0.7240, 0.4367))
sa <- cbind(sa, matrix(0,2,180))

# % reproduction
repro <- rbind(
    c(0, 2.5925, 4.4312, 5.1403, 5.5190, 5.6633, 5.6010, 5.4577, 5.2904, 5.0613, 4.6970, 
      3.3577, 1.5946, 1.0817, 0.9666, 0.8333, 0.4689, 0.0709, 
      0, 0, 0, 0),
    c(0, 3.1975, 5.4563, 6.2996, 6.7372, 6.9030, 6.8210, 6.6100, 6.1962, 5.1653, 4.1837, 
      3.6029, 3.1023, 2.4799, 1.6909, 1.1750, 1.0148, 0.9096, 
      0.7821, 0.6430, 0.5000, 0.3531)
    )
repro <- cbind(repro, matrix(0,2,178))


# % relative attack rate on the different instars
# % from Ives et al 1999
relatt_by_instar <- c(0.12, 0.27, 0.39, 0.16, 0.06)

# % putting in the attack rates for different development rates
relattemp <- matrix(0,2,totStage)
for (j in 1:2) {
    counter <- 0
    for (i in 1:5) {
        relattemp[j,(counter+1):(counter+stageInts[j,i])] <- relatt_by_instar[i] * 
            matrix(1,1,stageInts[j,i])
        counter <- counter + stageInts[j,i]
    }
}
relatt <- relattemp

sexratio <- 0.5


a <- 2.5
k <- 0.0005
kp <- 0.0006
kk <- 0.1811
h <- 0.0363

sw <- 0.55
s1 <- 0
s2 <- 0
s3 <- 0
rho <- 2 / (1 + exp(-s3)) - 1
mumDetect <- 0.3923
sm <- 33.5455

# % These are the survivals of singly attacked and multiply attacked
# % resistant aphids
fresist <- c(0.9, 0.6)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % set up Leslie matrices

# % resistant clones
L <- matrix(0,totStage,totStage)
juvTime <- sum(stageInts[clone[1,1],1:(end-1)])
# WHAT DOES -1 BELOW DO??
LL <- diag(c(sj[clone[1,1]] * matrix(1,1,juvTime), 
             sa[clone[1,2], 1:(totStage-juvTime-1)]), -1)
LL[1,(juvTime+1):(juvTime+stageInts[clone[1,1],end])] <- 
    repro[clone[1,2],1:(stageInts[clone[1,1],end])]
L <- LL

# Left off on line 87 of HighTunnelExptSim_14Aug17.m
