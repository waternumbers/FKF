#' Fast Kalman Smoother
#'
#' This function can be run after running \code{\link{fkf}} to produce
#' "smoothed" estimates of the state variable \code{a(t)} and it's variance
#' \code{V(t)}. Unlike the output of the filter, these estimates are conditional
#' on the entire data rather than only the past, that is it estimates \eqn{E[at
#' | y1,\ldots,yn]} and \eqn{V[at | y1,\ldots,yn]}.
#'
#' @param FKFobj  An S3-object of class "fkf", returned by \code{\link{fkf}}.
#'
#' @return A list with the following elements:
#'
#'   \code{ahatt}  A \eqn{m \times (n + 1)}{m * (n + 1)}-matrix containing the
#'   smoothed state variables, i.e. \eqn{ahat_t = E(\alpha_t | y_{n}}\cr
#'   \code{Vt}  A \eqn{m \times m \times (n + 1)}{m * m * (n + 1)}-array
#'   containing the variances of \code{ahatt}, i.e. \eqn{V_t = Var(\alpha_t |
#'   y_{n}}\cr
#'
#' @examples
#' ## <--------------------------------------------------------------------------->
#' ## Example 1: ARMA(2, 1) model estimation.
#' ## <--------------------------------------------------------------------------->
#' ## This example shows how to fit an ARMA(2, 1) model using this Kalman
#' ## filter implementation (see also stats' makeARIMA and KalmanRun).
#' n <- 1000
#'
#' ## Set the AR parameters
#' ar1 <- 0.6
#' ar2 <- 0.2
#' ma1 <- -0.2
#' sigma <- sqrt(0.2)
#'
#' ## Sample from an ARMA(2, 1) process
#' a <- arima.sim(model = list(ar = c(ar1, ar2), ma = ma1), n = n,
#'                innov = rnorm(n) * sigma)
#'
#' ## Create a state space representation out of the four ARMA parameters
#' arma21ss <- function(ar1, ar2, ma1, sigma) {
#'   Tt <- matrix(c(ar1, ar2, 1, 0), ncol = 2)
#'   Zt <- matrix(c(1, 0), ncol = 2)
#'   ct <- matrix(0)
#'   dt <- matrix(0, nrow = 2)
#'   GGt <- matrix(0)
#'   H <- matrix(c(1, ma1), nrow = 2) * sigma
#'   HHt <- tcrossprod(H)
#'   a0 <- c(0, 0)
#'   P0 <- matrix(1e6, nrow = 2, ncol = 2)
#'   return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
#'               HHt = HHt))
#' }
#'
#' ## The objective function passed to 'optim'
#' objective <- function(theta, yt) {
#'   sp <- arma21ss(theta["ar1"], theta["ar2"], theta["ma1"], theta["sigma"])
#'   ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
#'              Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
#'   return(-ans$logLik)
#' }
#'
#' theta <- c(ar = c(0, 0), ma1 = 0, sigma = 1)
#' fit <- optim(theta, objective, yt = rbind(a), hessian = TRUE)
#' fit
#'
#' ## Confidence intervals
#' rbind(fit$par - qnorm(0.975) * sqrt(diag(solve(fit$hessian))),
#'       fit$par + qnorm(0.975) * sqrt(diag(solve(fit$hessian))))
#'
#' ## Filter the series with estimated parameter values
#' sp <- arma21ss(fit$par["ar1"], fit$par["ar2"], fit$par["ma1"], fit$par["sigma"])
#' ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
#'            Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = rbind(a))
#' smooth <- fks(ans)
#'
#' ## Compare the filtered series with the realization
#' plot(ans, at.idx = NA, att.idx = 1, CI = NA)
#' lines(a, lty = "dotted")
#'
#' ## Compare the smoothed series with the realization
#' ##plot(smooth$ahatt[1,], col=1, type='l', lty=1)
#' ##lines(a, lty='dotted')
#'
#'
#' ## <--------------------------------------------------------------------------->
#' ## Example 2: Local level model for the Nile's annual flow.
#' ## <--------------------------------------------------------------------------->
#' ## Transition equation:
#' ## alpha[t+1] = alpha[t] + eta[t], eta[t] ~ N(0, HHt)
#' ## Measurement equation:
#' ## y[t] = alpha[t] + eps[t], eps[t] ~  N(0, GGt)
#'
#' y <- Nile
#' y[c(3, 10)] <- NA  # NA values can be handled
#'
#' ## Set constant parameters:
#' dt <- ct <- matrix(0)
#' Zt <- Tt <- matrix(1)
#' a0 <- y[1]            # Estimation of the first year flow
#' P0 <- matrix(100)     # Variance of 'a0'
#'
#' ## Estimate parameters:
#' fit.fkf <- optim(c(HHt = var(y, na.rm = TRUE) * .5,
#'                    GGt = var(y, na.rm = TRUE) * .5),
#'                  fn = function(par, ...)
#'                    -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
#'                  yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
#'                  Zt = Zt, Tt = Tt)
#'
#' ## Filter Nile data with estimated parameters:
#' fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]),
#'                GGt = matrix(fit.fkf$par[2]), yt = rbind(y))
#' ##fks.obj <- fks(fkf.obj)
#'
#' ## Compare with the stats' structural time series implementation:
#' fit.stats <- StructTS(y, type = "level")
#'
#' fit.fkf$par
#' fit.stats$coef
#'
#' ## Plot the flow data together with fitted local levels:
#' ##plot(y, main = "Nile flow")
#' ##lines(fitted(fit.stats), col = "green")
#' ##lines(ts(fkf.obj$att[1, ], start = start(y), frequency = frequency(y)), col = "blue")
#' ##lines(ts(fks.obj$ahatt[1,], start = start(y), frequency = frequency(y)), col = "red")
#' ##legend("top", c("Nile flow data", "Local level (StructTS)", "Local level (fkf)",
#' ##       "Local level (fks)"),
#' ##       col = c("black", "green", "blue", "red"), lty = 1)
#'
#' @export
fks <- function (FKFobj) {
  if (class(FKFobj) != 'fkf') stop('Input must be an object of class FKF')
  if (FKFobj$status[1] != 0 || FKFobj$status[2] != 0) {
    stop('Smoothing requires successful inversion of Ft for all t')
  }
  yt <- FKFobj$yt
  vt <- FKFobj$vt
  Ftinv <- FKFobj$Ftinv
  n <- dim(Ftinv)[3]
  Kt <- FKFobj$Kt
  at <- FKFobj$at[,1:n]
  Pt <- FKFobj$Pt[,,1:n]
  Zt <- FKFobj$Zt
  Tt <- FKFobj$Tt


  ans <- .Call("FKS", yt, Zt, vt, Tt, Kt, Ftinv, at, Pt, PACKAGE = "FKF")
  return(ans)
}
