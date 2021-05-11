#' Plotting fks objects
#'
#' Plotting method for objects of class \code{\link{fks}}. This function
#' provides tools visualisation of the state vector of the Kalman smoother output
#'
#' @param x The output of \code{\link{fks}}.
#' @param CI The confidence interval in case \code{type == "state"}. Set
#' \code{CI} to \code{NA} if no confidence interval shall be plotted.
#' @param ahatt.idx An vector giving the indexes of the predicted state variables
#' which shall be plotted if \code{type == "state"}.
#' @param \dots Arguments passed to either \code{\link{plot}},
#' \code{\link{qqnorm}}, \code{\link{qqplot}} or \code{\link{acf}}.
#'
#' @details
#' The state variables are plotted. By the argument \code{ahatt.idx}, the user can specify
#' which of the smoothed (\eqn{a_{t|n}}{a(t|n)}) state variables will be drawn.
#'
#' @return
#'
#' NULL
#' 
#' @seealso \code{\link{fks}}
#'
#' @examples
#' ## <--------------------------------------------------------------------------->
#' ## Example 3: Local level model for the treering data
#' ## <--------------------------------------------------------------------------->
#' ## Transition equation:
#' ## alpha[t+1] = alpha[t] + eta[t], eta[t] ~ N(0, HHt)          
#' ## Measurement equation:
#' ## y[t] = alpha[t] + eps[t], eps[t] ~  N(0, GGt)
#'
#' y <- treering
#' y[c(3, 10)] <- NA  # NA values can be handled
#'
#' ## Set constant parameters:
#' dt <- ct <- matrix(0)
#' Zt <- Tt <- array(1,c(1,1,1))
#' a0 <- y[1]            # Estimation of the first width
#' P0 <- matrix(100)     # Variance of 'a0'
#'
#' ## Estimate parameters:
#' fit.fkf <- optim(c(HHt = var(y, na.rm = TRUE) * .5,
#'                    GGt = var(y, na.rm = TRUE) * .5),
#'                  fn = function(par, ...)
#'                  -fkf(HHt = array(par[1],c(1,1,1)), GGt = array(par[2],c(1,1,1)), ...)$logLik,
#'                  yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
#'                  Zt = Zt, Tt = Tt)
#'
#' ## Filter tree ring data with estimated parameters:
#' fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = array(fit.fkf$par[1],c(1,1,1)),
#'                GGt = array(fit.fkf$par[2],c(1,1,1)), yt = rbind(y))
#'
#' fks.obj <- fks(fkf.obj)
#' plot(fks.obj)
#' lines(as.numeric(y),col="blue")
#' 
#'
#' @keywords plot
#' @export
plot.fks <- function(x,
                     CI = 0.95,
                     ahatt.idx = 1:nrow(x$ahatt), ...){
    
    d <- nrow(x$ahatt)
    n <- ncol(x$ahatt)
    
    xlim <- 1:n
    ylim <- range(x$ahatt[ahatt.idx, ], na.rm = TRUE)
    
    plot(xlim, ylim = ylim, type = "n", xlab = "Index",
             ylab = "State variables", ...)
        
    if(any(is.na(ahatt.idx))){
        ahatt.legend <- character(0)
        col.ahatt <- character(0)
        ahatt.idx <- numeric(0)
    }else{
        col.ahatt <- rainbow(length(ahatt.idx))
        for(i in 1:length(ahatt.idx)){
            lines(xlim, x$ahatt[ahatt.idx[i],],
                  col = col.ahatt[i], lty = "dashed")
            if(is.finite(CI)){
                lines(xlim, x$ahatt[ahatt.idx[i],] + qnorm(0.5 - CI / 2) * sqrt(x$Vt[ahatt.idx[i], ahatt.idx[i], ]),
                      col = col.ahatt[i], lty = "dotted")
                lines(xlim, x$ahatt[ahatt.idx[i],] + qnorm(0.5 + CI / 2) * sqrt(x$Vt[ahatt.idx[i], ahatt.idx[i], ]),
                      col = col.ahatt[i], lty = "dotted")
            }
        }
        ahatt.legend <- paste("ahatt[", ahatt.idx, ",]", sep = "")
    }

    legend("topleft", legend = ahatt.legend,
           lty = rep("dashed", length(ahatt.idx)),
           col = col.ahatt, bg = "white")

    if(is.finite(CI)){
      legend("topright", legend = ahatt.legend,
             lty = rep("dotted", length(ahatt.idx)),
             col = col.ahatt,
             title = paste("Confidence interval:", CI), bg = "white")
    }
}

