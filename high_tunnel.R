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
clone <- rbind(c(2, 1), c(2, 2))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % lab parameters

totStage <- 32

# % this sets the time scale in terms of the numbers of days per instar and
# % mummy development time
stageInts <- rbind(c(2, 2, 2, 2, 19), c(1, 1, 1, 2, 23))
mumInts <- cbind(7, 3)

# % juvenile survival
sj <- cbind(0.9745, 0.9849)

# % adult survival
sa <- rbind(c(1.0000, 0.9949, 0.9818, 0.9534, 0.8805, 0.8367, 0.8532, 0.8786, 0.8823, 
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
relatt_by_instar <- cbind(0.12, 0.27, 0.39, 0.16, 0.06)

# % putting in the attack rates for different development rates
relattemp <- matrix(0,2,totStage)
for (j in 1:2) {
    counter <- 0
    for (i in 1:5) {
        relattemp[j,(counter+1):(counter+stageInts[j,i])] <- relatt_by_instar[i] %*% 
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
fresist <- cbind(0.9, 0.6)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % set up Leslie matrices

# % resistant clones
L <- matrix(0,totStage,totStage)
juvTime <- sum(stageInts[clone[1,1],1:(ncol(stageInts)-1)])
LL <- diag_ml(c(sj[,clone[1,1]] * matrix(1,1,juvTime), 
                sa[clone[1,2], 1:(totStage-juvTime-1)]), -1)
LL[1,(juvTime+1):(juvTime+stageInts[clone[1,1],ncol(stageInts)])] <- 
    repro[clone[1,2],1:(stageInts[clone[1,1],ncol(stageInts)])]
L <- LL



L_eigen <- eigen(L)
SAD <- L_eigen$vectors
r <- matrix(L_eigen$values, ncol = 1)
rmax <- max(abs(r))
Rr <- rmax

SADdist <- SAD[, abs(r) == rmax]
SADdist <- SADdist / sum(SADdist)
Lr <- L
if (!any_complex(SADdist)) SADdist <- as.numeric(SADdist)
SADdistr <- matrix(SADdist, ncol = 1)

# % susceptible clones
L <- matrix(0,totStage,totStage)
juvTime <- sum(stageInts[clone[2,1],1:(ncol(stageInts)-1)])
LL <- diag_ml(c(sj[clone[2,1]] %*% matrix(1,1,juvTime), 
                sa[clone[2,2],1:(totStage-juvTime-1)]), -1)
LL[1,(juvTime+1):(juvTime+stageInts[clone[2,1],ncol(stageInts)])] <- 
    repro[clone[2,2], 1:stageInts[clone[2,1],ncol(stageInts)]]
L <- LL

L_eigen <- eigen(L)
SAD <- L_eigen$vectors
r <- L_eigen$values
rmax <- max(abs(r))
Rs <- rmax


SADdist <- SAD[, abs(r) == rmax]
SADdist <- SADdist / sum(SADdist)
if (!any_complex(SADdist)) SADdist <- as.numeric(SADdist)
Ls <- L
SADdists <- matrix(SADdist, ncol = 1)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % field parameters: This is set up to have different harvesting patterns
# % between nfield fields.

# % number of fields
nfield <- 2

# % kill rate at harvesting
kill <- 0.05

# % dispersal rates between fields for aphids, adult wasps, and predators
da <- 0.05
dw <- 1
pred <- 0.8

# % initial densities of aphids and parasitoids
initx <- 20
inity <- 1

# % time between harvests
cycles <- 20
cyclelength <- 30
fractresist <- 0.05

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % run program

xrinit <- fractresist * initx * SADdistr %*% matrix(1,1,nfield)
xsinit <- (1-fractresist) * initx * SADdists %*% matrix(1,1,nfield)

ny <- sum(mumInts) + 1
yr <- 0 * xrinit[1:(ny-1),]
yr <- inity * rbind(yr, c(1, 1))
ys <- 0 * xsinit[1:(ny-1),]
ys <- inity * rbind(ys, c(1, 1))


Tmax <- cyclelength * (1 + cycles)
harvesttimes <- rbind(c(cyclelength * (1:1), cyclelength * (2:cycles)),
                      c(cyclelength * (1:1), cyclelength * (2:(cycles-1)) - cyclelength/2,
                        cyclelength * cycles))

xs <- xsinit
xr <- xrinit










HighTunnelExptSimfunct <- function(
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
            
            
            mm <- rbind(a * relatt[clone[1,1],] * (yr[nrow(yr),i] + ys[nrow(ys),i]) / 
                            (h * sxyt + 1))
            AA <- (1 + t(mm) / kk)
            As <- AA^(-kk)
            
            mm <- rbind(a * relatt[clone[2,1],] * (yr[nrow(yr),i] + ys[nrow(ys),i]) / 
                            (h * sxyt + 1))
            AA <- (1 + t(mm) / kk)
            Ar <- AA^(-kk) + fresist[1] * t(mm) * AA^(-kk-1) + fresist[2] *
                (1-(AA^(-kk) + t(mm) * AA^(-kk-1)))
            
            xtr <- (pred * Kt * Ar) * as.matrix(LLr %*% xr[,i])
            xts <- (pred * Kt * As) * as.matrix(LLs %*% xs[,i])
            
            yt <- cbind(yr[,i])
            y <- cbind(yr[,i])
            x <- cbind(xr[,i])
            yt[length(yt)] <- sw * y[length(y)] + sexratio * y[(length(y)-1)]
            yt[(mumInts[1]+2):(length(yt)-1)] <- pred * y[(mumInts[1]+1):(length(y)-2)]
            yt[2:(mumInts[1]+1)] <- Kpt * sj[clone[1,1]] * y[1:mumInts[1]]
            yt[1] <- (Kpt * t(matrix(1,totStage,1) - Ar)) %*% as.matrix(LLr %*% x)
            ytr <- yt
            
            yt <- cbind(ys[,i])
            y <- cbind(ys[,i])
            x <- cbind(xs[,i])
            yt[length(yt)] <- sw * y[length(y)] + sexratio * y[(length(y)-1)]
            yt[(mumInts[1]+2):(length(yt)-1)] <- pred * y[(mumInts[1]+1):(length(y)-2)]
            yt[2:(mumInts[1]+1)] <- Kpt * sj[clone[2,1]] * y[1:mumInts[1]]
            yt[1] <- (Kpt * t(matrix(1,totStage,1) - As)) %*% as.matrix(LLs %*% x)
            yts <- yt
            
            if (t %in% harvesttimes[i,1:(ncol(harvesttimes)-1)]) {
                # Tony used a matrix for the end of a sequence (1:mumInts).
                # The behavior of this in matlab is to just end at the first element of 
                # the matrix (i.e., 1:mumInts[1,1])
                mumInts_ <- mumInts[1,1]
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
        
        dispersing <- da * cbind(rowMeans(xr[(sum(stageInts[clone[1,1],1:4])+1):nrow(xr),]))
        xr[(sum(stageInts[clone[1,1],1:4])+1):nrow(xr),] <- (1-da) *
            xr[(sum(stageInts[clone[1,1],1:4])+1):nrow(xr),] + 
            dispersing %*% matrix(1,1,nfield)
        
        dispersing <- da * cbind(rowMeans(xs[(sum(stageInts[clone[2,1],1:4])+1):nrow(xs),]))
        xs[(sum(stageInts[clone[2,1],1:4])+1):nrow(xs),] <- (1-da) *
            xs[(sum(stageInts[clone[2,1],1:4])+1):nrow(xs),] + 
            dispersing %*% matrix(1,1,nfield)
        
        dispersingw <- dw * mean(yr[nrow(yr),])
        yr[nrow(yr),] <- ((1-dw) * yr[nrow(yr),] + dispersingw)
        
        dispersingw <- dw * mean(ys[nrow(ys),])
        ys[nrow(ys),] <- ((1-dw) * ys[nrow(ys),] + dispersingw)
        
        nap <- nx + mumInts[1]
        
        yy <- rbind(xr, yr)
        Xr[t,] <- colSums(yy[1:nap,])
        Yr[t,] <- colSums(yy[(nx+1):nap,]) / colSums(yy[1:nap,])
        
        yy <- rbind(xs, ys)
        Xs[t,] <- colSums(yy[1:nap,])
        Ys[t,] <- colSums(yy[(nx+1):nap,]) / colSums(yy[1:nap,])
    }
    

    out_list <- list(Xr = Xr, Xs = Xs, Yr = Yr, Ys = Ys, xr = xr, xs = xs, 
                     yr = yr, ys = ys)
    
    return(out_list)
}



out_list <- HighTunnelExptSimfunct(xr,xs,yr,ys,a,fresist,k,kp,kk,h,sw,s1,s2,rho,sm,
                                   Tmax,nfield,kill,harvesttimes,da,dw,pred)
# Assigning out_list values to global ones for objects Xr, Xs, Yr, Ys, xr, xs, yr, & ys
invisible(
    lapply(names(out_list), 
           function(NAME) {
               eval(parse(text = paste0(NAME, ' <<- out_list$', NAME)))
           }))


Ymax <- 1.2 * max(max(Xr+Xs))
Ymax2 <- 1.2 * max(max(Xr+Xs))
range <- c(0, Tmax, 0, Ymax)
Tmin <- 1

# Figure 1
plot(1:(Tmax-Tmin+1),Xr[Tmin:Tmax,1]+Xs[Tmin:Tmax,1], type = 'l', col = 'dodgerblue')
lines(1:(Tmax-Tmin+1),Xr[Tmin:Tmax,2]+Xs[Tmin:Tmax,2], col = 'dodgerblue')
lines(1:(Tmax-Tmin+1),Ymax2 * (Yr[Tmin:Tmax,1]+Ys[Tmin:Tmax,1]), col = 'firebrick')
lines(1:(Tmax-Tmin+1),Ymax2 * (Yr[Tmin:Tmax,2]+Ys[Tmin:Tmax,2]), col = 'firebrick', lty = 2)
lines(1:(Tmax-Tmin+1),Ymax2*(Xr[Tmin:Tmax,1]+Xr[Tmin:Tmax,2])/
         (Xr[Tmin:Tmax,1]+Xr[Tmin:Tmax,2]+Xs[Tmin:Tmax,1]+Xs[Tmin:Tmax,2]), col = 'black')