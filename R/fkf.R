#' Fast Kalman filter
#'
#' This function allows for fast and flexible Kalman filtering. Both, the
#' measurement and transition equation may be multivariate and parameters
#' are allowed to be time-varying. In addition \dQuote{NA}-values in the
#' observations are supported. \code{fkf} wraps the \code{C}-function
#' \code{FKF} which fully relies on linear algebra subroutines contained
#' in BLAS and LAPACK.
#'
#' @section Usage:
#' \code{fkf(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt, check.input = TRUE)}
#'
#' @param a0 A \code{vector} giving the initial value/estimation of the state variable.
#' @param P0 A \code{matrix} giving the variance of \code{a0}.
#' @param dt A \code{matrix} giving the intercept of the transition equation (see \bold{Details}).
#' @param ct A \code{matrix} giving the intercept of the measurement equation (see \bold{Details}).
#' @param Tt An \code{array} giving the factor of the transition equation (see \bold{Details}).
#' @param Zt An \code{array} giving the factor of the measurement equation (see \bold{Details}).
#' @param HHt An \code{array} giving the variance of the innovations of the transition equation (see \bold{Details}).
#' @param GGt An \code{array} giving the variance of the disturbances of the measurement equation (see \bold{Details}).
#' @param yt A \code{matrix} containing the observations. \dQuote{NA}-values are allowed (see \bold{Details}).
#' @param check.input A \code{logical} stating whether the input shall be checked for consistency (\dQuote{storage.mode}, \dQuote{class}, and dimensionality, see \bold{Details}). This should always be 'TRUE' unless the performance becomes crucial and correctness of the arguments concerning dimensions, class and storage.mode is ensured.
#'
#' @details
#' \strong{State space form}
#'
#' The following notation is closest to the one of Koopman et al.
#' The state space model is represented by the transition equation and
#' the measurement equation. Let \eqn{m}{m} be the dimension of the state
#' variable, \eqn{d}{d} be the dimension of the observations, and \eqn{n}
#' the number of observations. The transition equation and the
#' measurement equation are given by
#' \deqn{\alpha_{t + 1} = d_t + T_t \cdot \alpha_t + H_t \cdot \eta_t}{alpha(t + 1) = d(t) + T(t) alpha(t) + H(t) * eta(t)} 
#' \deqn{y_t = c_t + Z_t \cdot \alpha_t + G_t \cdot \epsilon_t,}{y(t) = c(t) + Z(t) alpha(t) + G(t) * epsilon(t),}
#' where \eqn{\eta_t}{eta(t)} and \eqn{\epsilon_t}{epsilon(t)} are iid
#' \eqn{N(0, I_m)}{N(0, I(m))} and iid \eqn{N(0, I_d)}{N(0, I(d))},
#' respectively, and \eqn{\alpha_t}{alpha(t)} denotes the state
#' variable. The parameters admit the following dimensions:
#'
#' \tabular{lll}{
#' \eqn{a_t \in R^m} \tab \eqn{d_t \in R^m} \tab \eqn{eta_t \in R^m} \cr
#' \eqn{T_t \in R^{m \times m}}{d[t] \in R^(m * m)} \tab \eqn{H_t \in R^{m \times m}}{d[t] \in R^(m * m)} \tab \cr
#' \eqn{y_t \in R^d}{y[t] in R^d} \tab \eqn{c_t \in R^d}{c[t] \in R^d} \tab \eqn{\epsilon_t \in R^d}{epsilon[t] \in R^d} \cr
#' \eqn{Z_t \in R^{d \times m}}{Z[t] \in R^(d * m)} \tab \eqn{G_t \in R^{d \times d}}{G[t] \in R^(d * d)} \tab 
#' }
#'
#' Note that \code{fkf} takes as input \code{HHt} and \code{GGt} which
#' corresponds to \eqn{H_t H_t'}{H[t] \%*\% t(H[t])} and \eqn{G_t G_t'}{G[t] \%*\% t(G[t])}. 
#'
#' % <------------------------------------->
#' \strong{Iteration:}
#'
#' Let \code{i} be the loop variable. The filter iterations are
#' implemented the following way (in case of no NA's):
#'
#' Initialization:
#' \code{ if(i == 1)\{
#'   at[, i] = a0 
#'   Pt[,, i] = P0 
#' \} }
#'
#' Updating equations:\cr
#'   \code{vt[, i] = yt[, i] - ct[, i] - Zt[,,i] \%*\% at[, i]}\cr
#'   \code{Ft[,, i] = Zt[,, i] \%*\% Pt[,, i] \%*\% t(Zt[,, i]) + GGt[,, i]}\cr
#'   \code{Kt[,, i] = Pt[,, i] \%*\% t(Zt[,, i]) \%*\% solve(Ft[,, i])}\cr
#'   \code{att[, i] = at[, i] + Kt[,, i] \%*\% vt[, i]}\cr
#'   \code{Ptt[, i] = Pt[,, i] - Pt[,, i] \%*\% t(Zt[,, i]) \%*\% t(Kt[,, i])}
#'
#' Prediction equations:\cr
#'   \code{at[, i + 1] = dt[, i] + Tt[,, i] \%*\% att[, i]}\cr
#'   \code{Pt[,, i + 1] = Tt[,, i] \%*\% Ptt[,, i] \%*\% t(Tt[,, i]) + HHt[,, i]}
#'
#' Next iteration:\cr
#'   \code{i <- i + 1}\cr
#'   goto \dQuote{Updating equations}.
#'
#' % <------------------------------------->
#' \strong{NA-values:}
#'
#' NA-values in the observation matrix \code{yt} are supported.  If
#'   particular observations \code{yt[,i]} contain NAs, the NA-values are
#'   removed and the measurement equation is adjusted accordingly.  When
#'   the full vector \code{yt[,i]} is missing the Kalman filter reduces to
#'   a prediction step.
#'
#' % <------------------------------------->
#' \strong{Parameters:}
#'
#' The parameters can either be constant or deterministic
#'   time-varying. Assume the number of observations is \eqn{n}
#'   (i.e. \eqn{y = (y_t)_{t = 1, \ldots, n}, y_t = (y_{t1}, \ldots,
#'   y_{td})}{y = y[,1:n]}). Then, the parameters admit the following
#'   classes and dimensions:
#'
#' \tabular{ll}{
#'     \code{dt} \tab either a \eqn{m \times n}{m * n} (time-varying) or a \eqn{m \times 1}{m * 1} (constant) matrix. \cr
#'     \code{Tt} \tab either a \eqn{m \times m \times n}{m * m * n} or a \eqn{m \times m \times 1}{m * m * 1} array. \cr
#'     \code{HHt} \tab either a \eqn{m \times m \times n}{m * m * n} or a \eqn{m \times m \times 1}{m * m * 1} array. \cr
#'     \code{ct} \tab either a \eqn{d \times n}{d * n} or a \eqn{d \times 1}{d * 1} matrix. \cr
#'     \code{Zt} \tab either a \eqn{d \times m \times n}{d * m * n} or a \eqn{d \times m \times 1}{d * m * 1} array. \cr
#'     \code{GGt} \tab either a \eqn{d \times d \times n}{d * d * n} or a \eqn{d \times d \times 1}{d * d * 1} array. \cr
#'     \code{yt} \tab a \eqn{d \times n}{d * n} matrix.
#'   }
#'
#' If \code{check.input} is \code{TRUE} each argument will be checked for
#'   correctness of the dimensionality, storage mode, and
#'   class. \code{check.input} should always be \code{TRUE} unless the
#'   performance becomes crucial and correctness of the arguments
#'   concerning dimensions, class, and storage.mode is ensured.\cr
#'   Note that the class of the arguments if of importance. For instance,
#'   to check whether a parameter is constant the \code{dim} attribute is
#'   accessed. If, e.g., \code{Zt} is a constant, it could be a \eqn{d
#'   \times d}{d * d}-matrix. But the third dimension
#'   (i.e. \code{dim(Zt)[3]}) is needed to check for constancy. This
#'   requires \code{Zt} to be an \eqn{d \times d \times 1}{d * d *
#'   1}-array.
#'
#' % <------------------------------------->
#'   \strong{BLAS and LAPACK routines used:}
#'
#' The \R function \code{fkf} basically wraps the \code{C}-function
#'   \code{FKF}, which entirely relies on linear algebra subroutines
#'   provided by BLAS and LAPACK. The following functions are used:
#'
#' \tabular{rl}{
#'     BLAS: \tab \code{dcopy}, \code{dgemm}, \code{daxpy}. \cr
#'     LAPACK: \tab \code{dpotri}, \code{dpotrf}.
#'   }
#'
#' \code{FKF} is called through the \code{.Call} interface.  Internally,
#'   \code{FKF} extracts the dimensions, allocates memory, and initializes
#'   the \R-objects to be returned. \code{FKF} subsequently calls
#'   \code{cfkf} which performs the Kalman filtering.
#'
#' The only critical part is to compute the inverse of \eqn{F_t}{F[,,t]}
#'   and the determinant of \eqn{F_t}{F[,,t]}. If the inverse can not be
#'   computed, the filter stops and returns the corresponding message in
#'   \code{status} (see \bold{Value}). If the computation of the
#'   determinant fails, the filter will continue, but the log-likelihood
#'   (element \code{logLik}) will be \dQuote{NA}.
#'
#'   The inverse is computed in two steps:
#' First, the Cholesky factorization of \eqn{F_t}{F[,,t]} is
#'   calculated by \code{dpotrf}. Second, \code{dpotri} calculates the
#'   inverse based on the output of \code{dpotrf}. \cr
#'   The determinant of \eqn{F_t}{F[,,t]} is computed using again the
#'   Cholesky decomposition.
#' 
#' @return
#' An S3-object of class \dQuote{fkf}, which is a list with the following elements:
#' 
#' \tabular{rl}{
#'     \code{att} \tab A \eqn{m \times n}{m * n}-matrix containing the filtered state variables, i.e. \eqn{a_{t|t} = E(\alpha_t | y_t)}{att[,t] = E(alpha[t] | y[,t])}.\cr
#'     \code{at} \tab A \eqn{m \times (n + 1)}{m * (n + 1)}-matrix containing the predicted state variables, i.e. \eqn{a_t = E(\alpha_t | y_{t - 1})}{at[,t] = E(alpha[t] | y[,t - 1])}.\cr
#'     \code{Ptt} \tab A \eqn{m \times m \times n}{m * m * n}-array containing the variance of \code{att}, i.e. \eqn{P_{t|t} =  var(\alpha_t | y_t)}{Ptt[,,t] = var(alpha[t] | y[,t])}.\cr
#'     \code{Pt} \tab A \eqn{m \times m \times (n + 1)}{m * m * (n + 1)}-array containing the variances of \code{at}, i.e. \eqn{P_t = var(\alpha_t | y_{t - 1})}{Pt[,,t] = var(alpha[t] | y[,t - 1])}.\cr
#'     \code{vt} \tab A \eqn{d \times n}{d * n}-matrix of the prediction errors given by \eqn{v_t = y_t - c_t - Z_t a_t}{vt[,t] = yt[,t] - ct[,t] - Zt[,,t] \%*\% at[,t]}.\cr
#'     \code{Ft} \tab A \eqn{d \times d \times n}{d * d * n}-array which contains the variances of \code{vt}, i.e. \eqn{F_t = var(v_t)}{Ft[,,t] = var(v[,t])}.\cr
#'     \code{Kt} \tab A \eqn{m \times d \times n}{m * d * n}-array containing the \dQuote{Kalman gain} (ambiguity, see calculation above). \cr
#'     \code{logLik} \tab The log-likelihood. \cr
#'     \code{status} \tab A vector which contains the status of LAPACK's \code{dpotri} and \code{dpotrf}. \eqn{(0, 0)} means successful exit.\cr
#'   \code{sys.time} \tab The time elapsed as an object of class \dQuote{proc_time}.
#' }
#' 
#' The first element of both \code{at} and \code{Pt} is filled with the
#' function arguments \code{a0} and \code{P0}, and the last, i.e. the (n +
#' 1)-th, element of \code{at} and \code{Pt} contains the predictions \cr
#' \eqn{at[,n + 1] = E(\alpha_{n + 1} | y_n)}{at[,n + 1] = E(alpha[n + 1] | y[,n])} and \cr
#' \eqn{Pt[,,n + 1] = var(\alpha_{n + 1} | y_n)}{Pt[,,n + 1] = var(alpha[n + 1] | y[,n])}.
#'
#' @section References:
#'   Harvey, Andrew C. (1990). \emph{Forecasting, Structural Time Series
#'   Models and the Kalman Filter}.  Cambridge University Press.
#'
#' Hamilton, James D. (1994). \emph{Time Series Analysis}.  Princeton
#' University Press.
#'
#' Koopman, S. J., Shephard, N., Doornik, J. A. (1999).
#' \emph{Statistical algorithms for models in state space using SsfPack
#' 2.2}. Econometrics Journal, Royal Economic Society, vol. 2(1), pages
#' 107-160.
#'
#' @seealso \code{\link[=plot.fkf]{plot}} to visualize and analyze \code{fkf}-objects, \code{\link{KalmanRun}} from the stats package, function \code{dlmFilter} from package \code{dlm}.
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
#'     Tt <- matrix(c(ar1, ar2, 1, 0), ncol = 2)
#'     Zt <- matrix(c(1, 0), ncol = 2)
#'     ct <- matrix(0)
#'     dt <- matrix(0, nrow = 2)
#'     GGt <- matrix(0)
#'     H <- matrix(c(1, ma1), nrow = 2) * sigma
#'     HHt <- H %*% t(H)
#'     a0 <- c(0, 0)
#'     P0 <- matrix(1e6, nrow = 2, ncol = 2)
#'     return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
#'                 HHt = HHt))
#' }
#'
#' ## The objective function passed to 'optim'
#' objective <- function(theta, yt) {
#'     sp <- arma21ss(theta["ar1"], theta["ar2"], theta["ma1"], theta["sigma"])
#'     ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
#'                Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
#'     return(-ans$logLik)
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
#'
#' ## Compare the prediction with the realization
#' plot(ans, at.idx = 1, att.idx = NA, CI = NA)
#' lines(a, lty = "dotted")
#'
#' ## Compare the filtered series with the realization
#' plot(ans, at.idx = NA, att.idx = 1, CI = NA)
#' lines(a, lty = "dotted")
#'
#' ## Check whether the residuals are Gaussian
#' plot(ans, type = "resid.qq")
#'
#' ## Check for linear serial dependence through 'acf'
#' plot(ans, type = "acf")
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
#'                  -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
#'                  yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
#'                  Zt = Zt, Tt = Tt)
#'
#' ## Filter Nile data with estimated parameters:
#' fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]),
#'                GGt = matrix(fit.fkf$par[2]), yt = rbind(y))
#'
#' ## Compare with the stats' structural time series implementation:
#' fit.stats <- StructTS(y, type = "level")
#'
#' fit.fkf$par
#' fit.stats$coef
#'
#' ## Plot the flow data together with fitted local levels:
#' plot(y, main = "Nile flow")
#' lines(fitted(fit.stats), col = "green")
#' lines(ts(fkf.obj$att[1, ], start = start(y), frequency = frequency(y)), col = "blue")
#' legend("top", c("Nile flow data", "Local level (StructTS)", "Local level (fkf)"),
#'        col = c("black", "green", "blue"), lty = 1)
#'
#' @keywords algebra models multivariate
#' @export
fkf <- function(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt, check.input = TRUE) {
    if(any(c(missing(a0), missing(P0), missing(dt), missing(ct), missing(Tt),
             missing(Zt), missing(HHt), missing(GGt), missing(yt)))){
        
        stop("None of the input arguments 'a0', 'P0', 'dt', 'ct', 'Tt', 'Zt',",
             "'HHt', 'GGt', and 'yt' must be missing.")
    }

    ## 'check.input' should always be 'TRUE' unless the performance
    ## becomes crucial and correctness of the arguments concerning
    ## dimensions, class and storage.mode is ensured.
    if(check.input){

        ## Check the storage mode: Must be 'double' for all arguments
        stor.mode <- c(storage.mode(a0), storage.mode(P0), storage.mode(dt),
                       storage.mode(ct), storage.mode(Tt), storage.mode(Zt),
                       storage.mode(HHt), storage.mode(GGt), storage.mode(yt))

        names(stor.mode) <- c("a0", "P0", "dt", "ct", "Tt", "Zt",
                              "HHt", "GGt", "yt")

        if(any(stor.mode != "double")){
            stop("storage mode of variable(s) '",
                 paste(names(stor.mode)[stor.mode != "double"],
                       collapse = "', '"),
                 "' is not 'double'!\n", sep = "")
        }

        ## Check classes of arguments
        error.string <- ""
        if(!is.vector(a0)){
            error.string <- paste(error.string,
                                  "'a0' must be a vector!\n", sep = "")
        }

        if(!is.matrix(P0)){
            error.string <- paste(error.string,
                                  "'P0' must be of class 'matrix'!\n", sep = "")
        }

        if(!is.matrix(dt) && !is.vector(dt)){
            error.string <- paste(error.string,
                                  "'dt' must be of class 'vector' or 'matrix'!\n", sep = "")
        }else if(is.vector(dt)){
            dt <- as.matrix(dt)
        }

        if(!is.matrix(ct) && !is.vector(ct)){
            error.string <- paste(error.string,
                                  "'ct' must be of class 'vector' or 'matrix'!\n", sep = "")
        }else if(is.vector(ct)){
            ct <- as.matrix(ct)
        }

        if(!is.array(Tt)){
            error.string <- paste(error.string,
                                  "'Tt' must be of class 'matrix' or 'array'!\n", sep = "")
        }else if(is.matrix(Tt)){
            Tt <- array(Tt, c(dim(Tt), 1))
        }

        if(!is.array(Zt)){
            error.string <- paste(error.string,
                                  "'Zt' must be of class 'matrix' or 'array'!\n", sep = "")
        }else if(is.matrix(Zt)){
            Zt <- array(Zt, c(dim(Zt), 1))
        }

        if(!is.array(HHt)){
            error.string <- paste(error.string,
                                  "'HHt' must be of class 'matrix' or 'array'!\n", sep = "")
        }else if(is.matrix(HHt)){
            HHt <- array(HHt, c(dim(HHt), 1))
        }

        if(!is.array(GGt)){
            error.string <- paste(error.string,
                                  "'GGt' must be of class 'matrix' or 'array'!\n", sep = "")
        }else if(is.matrix(GGt)){
            GGt <- array(GGt, c(dim(GGt), 1))
        }

        if(!is.matrix(yt)){
            error.string <- paste(error.string,
                                  "'yt' must be of class 'matrix'!\n", sep = "")
        }
        if(error.string != ""){
            stop(error.string)
        }
        ## Check compatibility of dimensions
        n <- ncol(yt)
        d <- nrow(yt)
        m <- length(a0)

        if(dim(P0)[2] != m | dim(P0)[1] != m | dim(dt)[1] != m |
           dim(Zt)[2] != m | dim(HHt)[1] != m | dim(HHt)[2] != m  |
           dim(Tt)[1] != m  | dim(Tt)[2] != m)
        {
            stop("Some of dim(P0)[2], dim(P0)[1], dim(dt)[1],\n",
                 "dim(Zt)[2], dim(HHt)[1], dim(HHt)[2],\n",
                 "dim(Tt)[1] or dim(Tt)[2] is/are not equal to 'm'!\n")
        }

        if((dim(dt)[2] != n && dim(dt)[2] != 1) |
           (dim(ct)[2] != n && dim(ct)[2] != 1) |
           (dim(Tt)[3] != n && dim(Tt)[3] != 1) |
           (dim(Zt)[3] != n && dim(Zt)[3] != 1) |
           (dim(HHt)[3] != n && dim(HHt)[3] != 1) |
           (dim(GGt)[3] != n && dim(GGt)[3] != 1) |
           dim(yt)[2]  != n)
        {
            stop("Some of dim(dt)[2], dim(ct)[2], dim(Tt)[3],\n",
                 "dim(Zt)[3], dim(HHt)[3], dim(GGt)[3] or\n",
                 "dim(yt)[2] is/are neither equal to 1 nor equal to 'n'!\n")
        }

        if(dim(ct)[1] != d | dim(Zt)[1] != d  |
           dim(GGt)[1]!= d  | dim(GGt)[2] != d  | dim(yt)[1] != d)
        {
            stop("Some of dim(ct)[1], dim(Zt)[1], dim(GGt)[1],\n",
                 "dim(GGt)[2] or dim(yt)[1] is/are not equal to 'd'!\n")

        }

    }

    time.0 <- proc.time()

    ans <-.Call("FKF", a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt, PACKAGE = "FKF")

    ans$sys.time <- proc.time() - time.0

    return(ans)
}

