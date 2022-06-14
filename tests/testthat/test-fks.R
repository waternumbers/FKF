test_that("FKS matches simple implimentation", {
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
    Y <- Z + matrix(rnorm(2*n,0,sqrt(2)),2,n) ## add noise
    
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

    
    ## filter estimates using fkf
    fkf.obj <- fkf(a0,P0,dt,ct,Tt,Zt,HHt,GGt,Y)

    ## simple smoother implimentation - direct coding of the documentation
    d <- nrow(fkf.obj$vt)
    n <- ncol(fkf.obj$vt)
    r <- matrix(0,d,1)
    N <- matrix(0,d,d)

    ## implimented as for Durbinand Koopman paper
    Vt <- fkf.obj$Pt[,,1:n]
    mt <- fkf.obj$at[,1:n]
    for(ii in (n:1)){
        
        ## change r and N
        idx <- is.finite(fkf.obj$vt[,ii])
        K <- fkf.obj$Kt[,,ii]; K <- K[,idx,drop=FALSE]
        ZZ <- fkf.obj$Zt[,,1]; ZZ <- ZZ[idx,,drop=FALSE]
        L <- as.matrix(fkf.obj$Tt[,,1]) %*% (diag(d) - K%*%ZZ) #fkf.obj$Kt[,idx,ii]%*%fkf.obj$Zt[idx,,1]
        r <- t(ZZ) %*% fkf.obj$Ftinv[idx,idx,ii] %*% matrix(fkf.obj$vt[idx,ii]) + t(L)%*%r
        N <- t(ZZ)%*% fkf.obj$Ftinv[idx,idx,ii] %*% ZZ + t(L)%*%N%*%L
        
        ## apply - recall at start Vt[,,ii] is Pt[,,ii] and mt[,ii] are filtered estimates filtered estimate
        mt[,ii] <- mt[,ii] + Vt[,,ii]%*%r
        Vt[,,ii] <- Vt[,,ii] - Vt[,,ii]%*%N%*%Vt[,,ii]
        ##Nrec[,,ii] <- N
        

    }
    

    ## run fks
    fks.obj <- fks(fkf.obj)
    
    ## tests
    testthat::expect_equal( mt,fks.obj$ahatt )
    testthat::expect_equal( Vt[,,1:1000], fks.obj$Vt )

})

## Multiple missing inputs
test_that("FKS handles multiple missing data in transpsed matrix", {

    fk <- readRDS(testthat::test_path("data", "multi_missing_obs.rds"))
    out <- fks(fk)
    expect_s3_class(out, "fks")
})
