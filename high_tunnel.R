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
            
            
            # mm=a*relatt(clone(1,1),:)*(yr(end,i)+ys(end,i))/(h*sxyt+1);
            # AA=(1+mm'/kk);
            # As=AA.^(-kk);
                
            # mm=a*relatt(clone(2,1),:)*(yr(end,i)+ys(end,i))/(h*sxyt+1);
            # AA=(1+mm'/kk);
            # Ar=AA.^(-kk) + fresist(1)*mm'.*AA.^(-kk-1) + fresist(2)*(1-(AA.^(-kk)+mm'.*AA.^(-kk-1)));

            
            
        }
    }
    
    
    
}

