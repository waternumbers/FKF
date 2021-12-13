test_that("FKF matches stats::StructTS", {
    ## this is the original test from RUnit
    ## Fit a l
    
    data <- as.numeric(treering)
    data[c(3, 10:31)] <- NA
    
    fit.stats <- StructTS(data, type = "level")
    
    yt <- rbind(data)
    a0 <- yt[1]
    P0 <- matrix(100)
    dt <- ct <- matrix(0)
    Zt <- Tt <- matrix(1)
    HHt <- var(t(yt), na.rm = TRUE)
    GGt <- var(t(yt), na.rm = TRUE) * 0.1
    
    StructTS.level <- function(par, ...){
        return(-fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik)
    }
    fit.fkf <- optim(c(HHt, GGt), StructTS.level, yt = yt,
                     a0 = a0, P0 = P0, dt = dt, ct = ct,
                     Zt = Zt, Tt = Tt)
    
    expect_equal(as.numeric(fit.fkf$par), as.numeric(fit.stats$coef), tolerance=0.01)
    ## from original RUnit code
    ##checkTrue(all((fit.fkf$par / fit.stats$coef - 1) < 0.01),
    ##          "Difference between 'FKF' and 'stats' implementation is greater than 1 percent!")
})

test_that("FKF matches simple implimentation", {
    ## This tests the code for handling of missing data
    ## Comparision to a simple implimentation in R
    
    ## create data from an integrated random walk
    n <- 1000
    Z <- matrix(NA,2,n)
    Z[1,1] <- rnorm(1) # initial condition of integral
    Z[2,] <- rnorm(n) # perturbations of random walk
    for(ii in 2:n){
        Z[2,ii] <- Z[2,ii] + Z[2,ii-1] # evolve random walk
        Z[1,ii] <- Z[1,ii-1] + Z[2,ii-1] # evolve integral
    }
    Y <- Z + matrix(rnorm(2*n,0,sqrt(0.2)),2,n) ## add noise
    
    Y[,20:50] <- NA
    Y[1, c(100:105,300,544)] <- NA
    Y[2, c(400:405,997,n)] <- NA
    
    ## Set up model - this is the true model
    dt <- ct <- array(0,c(2,1))
    Zt <- array(c(1,0,1,0),c(2,2,1))
    Tt <- array(c(1,0,1,1),c(2,2,1))
    a0 <- Y[,1]            # Estimation of the first width
    P0 <- diag(2)*100     # Variance of 'a0'
    HHt <- array(c(0,0,0,1),c(2,2,1))
    GGt <- array(c(0.2,0,0,0.2),c(2,2,1))
    
    ## simple filter - direct coding of documentation
    at <- array(NA,c(2,n+1))
    att <- array(NA,c(2,n))
    Pt <- array(NA,c(2,2,n+1))
    Ptt <- array(NA,c(2,2,n))
    at[,1] <- a0
    Pt[,,1] <- P0
    for(ii in 1:n){
        ## correct
        idx <- is.finite(Y[,ii])
        if(any(idx)){
            v = Y[idx, ii] - ct[idx,1] - Zt[idx,,1] %*% at[, ii]
            Z <- Zt[,,1]
            Z <- Z[idx,,drop=FALSE]
            Ft = Z %*% Pt[,, ii] %*% t(Z) + GGt[idx,idx, 1]
            Kt = Pt[,, ii] %*% t(Z) %*% solve(Ft)
            att[, ii] = at[, ii] + Kt %*% v
            Ptt[,, ii] = Pt[,, ii] - Pt[,, ii] %*% t(Z) %*% t(Kt)
        }else{
            att[, ii] = at[, ii]
            Ptt[,, ii] = Pt[,, ii]
        }
        
        ## predict
        at[,ii+1] <- dt[,1] + Tt[,,1]%*% att[,ii]
        Pt[,, ii + 1] = Tt[,, 1] %*% Ptt[,, ii] %*% t(Tt[,, 1]) + HHt[,, 1]
    }
    ## fkf function
    fkf.obj <- fkf(a0,P0,dt,ct,Tt,Zt,HHt,GGt,Y)

    ## tests
    testthat::expect_equal( at,fkf.obj$at )
    testthat::expect_equal( Pt, fkf.obj$Pt )
    testthat::expect_equal( att,fkf.obj$att )
    testthat::expect_equal( Ptt, fkf.obj$Ptt )

})
