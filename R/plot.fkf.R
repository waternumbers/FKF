#' Plotting fkf objects
#'
#' Plotting method for objects of class \code{\link{fkf}}. This function
#' provides tools for graphical analysis of the Kalman filter output:
#' Visualization of the state vector, QQ-plot of the individual
#' residuals, QQ-plot of the Mahalanobis distance, auto- as well as
#' crosscorrelation function of the residuals.
#'
#' @section usage:
#' \code{plot(x, type = c("state", "resid.qq", "qqchisq", "acf"),
#' CI = 0.95, at.idx = 1:nrow(x$at), att.idx = 1:nrow(x$att), \dots)}
#'
#' @param x The output of \code{\link{fkf}}.
#' @param type A string stating what shall be plotted (see \bold{Details}).
#' @param CI The confidence interval in case \code{type == "state"}. Set
#' \code{CI} to \code{NA} if no confidence interval shall be plotted.
#' @param at.idx An vector giving the indexes of the predicted state variables
#' which shall be plotted if \code{type == "state"}.
#' @param att.idx An vector giving the indexes of the filtered state variables
#' which shall be plotted if \code{type == "state"}.
#' @param \dots Arguments passed to either \code{\link{plot}},
#' \code{\link{qqnorm}}, \code{\link{qqplot}} or \code{\link{acf}}.
#'
#' @details
#' The argument \code{type} states what shall be plotted. \code{type}
#' must partially match one of the following:
#' \describe{
#'     \item{\code{state}}{The state variables are plotted. By the
#'       arguments \code{at.idx} and \code{att.idx}, the user can specify
#'       which of the predicted (\eqn{a_{t}}{at}) and filtered
#'       (\eqn{a_{t|t}}{att}) state variables will be drawn.}
#'       \item{\code{resid.qq}}{Draws a QQ-plot for each residual-series in\code{vt}.}
#'       \item{\code{qqchisq}}{A Chi-Squared QQ-plot will be drawn to graphically
#' 	test for multivariate normality of the residuals based on the
#' 	Mahalanobis distance.}
#'       \item{\code{acf}}{Creates a pairs plot with the autocorrelation
#' 	function (\code{\link{acf}}) on the diagonal panels and the
#' 	crosscorrelation function (\code{\link{ccf}}) of the residuals on the
#' 	off-diagnoal panels.}
#'     }
#'
#' @return
#'   Invisibly returns an list with components:
#'    \tabular{rl}{
#'   \code{distance} \tab The Mahalanobis distance of the residuals as a
#'   vector of length \eqn{n}. \cr
#'   \code{std.resid} \tab The standardized residuals as an \eqn{d \times
#'     n}{d * n}-matrix. It should hold that \eqn{std.resid_{ij} \; iid \sim N_d(0, I)}{std.resid[i,j] iid N_d(0,
#'       I)},
#' }
#' where \eqn{d} denotes the dimension of the data and \eqn{n} the number
#' of observations.
#'
#' @seealso \code{\link{fkf}}
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
#'                  Zt = Zt, Tt = Tt, check.input = FALSE)
#'
#' ## Filter tree ring data with estimated parameters:
#' fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = array(fit.fkf$par[1],c(1,1,1)),
#'                GGt = array(fit.fkf$par[2],c(1,1,1)), yt = rbind(y))
#'
#' ## Plot the width together with fitted local levels:
#' plot(y, main = "Treering data")
#' lines(ts(fkf.obj$att[1, ], start = start(y), frequency = frequency(y)), col = "blue")
#' legend("top", c("Treering data", "Local level"), col = c("black", "blue"), lty = 1)
#'
#' ## Check the residuals for normality:
#' plot(fkf.obj, type = "resid.qq")
#'
#' ## Test for autocorrelation:
#' plot(fkf.obj, type = "acf", na.action = na.pass)
#'
#' @keywords hplot
#' @export
plot.fkf <- function(x, type = c("state", "resid.qq", "qqchisq", "acf"),
                     CI = 0.95,
                     at.idx = 1:nrow(x$at), att.idx = 1:nrow(x$att), ...){
    type <- match.arg(type)
    
    d <- nrow(x$vt)
    n <- ncol(x$att)
    
    distance <- sapply(1:n, function(i, Ft, vt){
        if(any(is.na(vt[, i])))
            return(matrix(NA))
        else
            return(t(vt[,i]) %*% solve(Ft[,,i]) %*% vt[,i])
    },
    Ft = x$Ft, vt = x$vt)
    
    std.resid <- sapply(1:n, function(i, Ft, vt)
    {
        if(any(is.na(vt[, i])))
            return(matrix(NA, ncol = 1, nrow = nrow(vt)))
        else
            return(solve(t(chol(Ft[,,i]))) %*% vt[,i])
    },
    Ft = x$Ft, vt = x$vt)
    
    dim(std.resid) <- c(d, n) # In case d == 1, make std.resid to be a matrix.
    ##< ------------------------------------------------------------ >
    if(type == "state"){
        xlim <- 1:(n + 1)
        xlim.att <- 1:n                     # Nbr. of predictions - 1
        
        ylim <- range(x$att[att.idx, ], x$at[at.idx, ], na.rm = TRUE)
        plot(xlim, ylim = ylim, type = "n", xlab = "Index",
             ylab = "State variables", ...)
        
    if(any(is.na(at.idx))){
      at.legend <- character(0)
      col.at <- character(0)
      at.idx <- numeric(0)
    }else{
      col.at <- rainbow(length(at.idx))
      for(i in 1:length(at.idx)){
        lines(xlim, x$at[at.idx[i],],
              col = col.at[i], lty = "dashed")
        if(is.finite(CI)){
          lines(xlim, x$at[at.idx[i],] + qnorm(0.5 - CI / 2) * sqrt(x$Pt[at.idx[i], at.idx[i], ]),
                col = col.at[i], lty = "dotted")
          lines(xlim, x$at[at.idx[i],] + qnorm(0.5 + CI / 2) * sqrt(x$Pt[at.idx[i], at.idx[i], ]),
                col = col.at[i], lty = "dotted")
        }
      }
      at.legend <- paste("at[", at.idx, ",]", sep = "")
    }

    if(any(is.na(att.idx))){
      att.legend <- character(0)
      col.att <- character(0)
      att.idx <- numeric(0)
    }else{
      col.att <- rainbow(length(att.idx))
      for(i in 1:length(att.idx)){
        lines(xlim.att, x$att[att.idx[i],],
              col = col.att[i], lty = "solid")
        if(is.finite(CI)){
          lines(xlim.att, x$att[att.idx[i],] + qnorm(0.5 - CI / 2) * sqrt(x$Ptt[att.idx[i], att.idx[i], ]),
                col = col.att[i], lty = "dotdash")
          lines(xlim.att, x$att[att.idx[i],] + qnorm(0.5 + CI / 2) * sqrt(x$Ptt[att.idx[i], att.idx[i], ]),
                col = col.att[i], lty = "dotdash")
        }

      }
      att.legend <- paste("att[", att.idx, ",]", sep = "")
    }

    legend("topleft", legend = c(at.legend, att.legend),
           lty = c(rep("dashed", length(at.idx)),
             rep("solid", length(att.idx))),
           col = c(col.at, col.att), bg = "white")

    if(is.finite(CI)){
      legend("topright", legend = c(att.legend, at.legend),
             lty = c(rep("dotted", length(at.idx)),
               rep("dotdash", length(att.idx))),
             col = c(col.at, col.att),
             title = paste("Confidence interval:", CI), bg = "white")
    }
    ## < ------------------------------------------------------------ >
  }else if(type == "resid.qq"){
    if(d < 3){
      par(mfrow = c(d, 1))
    }else{
      par(mfrow = c(d %/% 3 + (d %% 3 > 0), 3))
    }
    for(i in 1:d){
      qqnorm(std.resid[i, ],
             main = paste("Normal Q-Q Plot of residuals 'vt[",i,",]'", sep = ""),
             ...)
    }
    ## < ------------------------------------------------------------ >
  }else if(type == "qqchisq"){
    q.distance <- qchisq(ppoints(n), df = d)[order(order(distance))]
    qqplot(q.distance, distance,
           main = expression(paste(chi[2]^n,
               " Q-Q Plot of transformed residuals Dt = t(vt) * Ft^(-1) * vt")),
           xlab = "Theoretical Quantiles",
           ylab = "Sample Quantiles", ...)
    ## < ------------------------------------------------------------ >
  }else if(type == "acf"){
    par(mfrow = c(d, d))
    for(i in 1:d){
      for(j in 1:d){
        if(i == j){
          acf(x$vt[i, ],
              main = paste("ACF: vt[", i,", ] & vt[", i,", ]", sep = ""), ...)
        }else{
          ccf(x$vt[i, ], x$vt[j, ],
              main = paste("CCF: vt[", i,", ] & vt[", j,", ]", sep = ""), ...)
        }
      }
    }

  }
  invisible(list(std.resid = std.resid, distance = distance))
}

