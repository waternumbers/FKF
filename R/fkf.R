#' Fast Kalman filter
#'
#' @description This function allows for fast and flexible Kalman filtering. Both, the
#' measurement and transition equation may be multivariate and parameters
#' are allowed to be time-varying. In addition \dQuote{NA}-values in the
#' observations are supported. \code{fkf} wraps the \code{C}-function
#' \code{FKF} which fully relies on linear algebra subroutines contained
#' in BLAS and LAPACK.
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
# #' @param check.input A \code{logical} stating whether the input shall be checked for consistency (\dQuote{storage.mode}, \dQuote{class}, and dimensionality, see \bold{Details}). This input is depreciated and will be removed in a future version, checks are always made.
#'
#' @return
#' An S3-object of class \dQuote{fkf}, which is a list with the following elements:
#' 
#' \tabular{rl}{
#'     \code{att} \tab A \eqn{m \times n}{m * n}-matrix containing the filtered state variables, i.e. att[,t] = \eqn{a_{t|t}}{a(t|t)}.\cr
#'     \code{at} \tab A \eqn{m \times (n + 1)}{m * (n + 1)}-matrix containing the predicted state variables, i.e. at[,t] = \eqn{a_t}{a(t)}.\cr
#'     \code{Ptt} \tab A \eqn{m \times m \times n}{m * m * n}-array containing the variance of \code{att}, i.e. Ptt[,,t] = \eqn{P_{t|t}}{P(t|t)}.\cr
#'     \code{Pt} \tab A \eqn{m \times m \times (n + 1)}{m * m * (n + 1)}-array containing the variances of \code{at}, i.e. Pt[,,t] = \eqn{P_t}{P(t)}.\cr
#'     \code{vt} \tab A \eqn{d \times n}{d * n}-matrix of the prediction errors i.e. vt[,t] = \eqn{v_t}{v(t)}.\cr
#'     \code{Ft} \tab A \eqn{d \times d \times n}{d * d * n}-array which contains the variances of \code{vt}, i.e. Ft[,,t] = \eqn{F_t}{F(t)}.\cr
#'     \code{Kt} \tab A \eqn{m \times d \times n}{m * d * n}-array containing the \dQuote{Kalman gain} i.e. Kt[,,t] = \eqn{k_t}{K(t)}.\cr
#'     \code{logLik} \tab The log-likelihood. \cr
#'     \code{status} \tab A vector which contains the status of LAPACK's \code{dpotri} and \code{dpotrf}. \eqn{(0, 0)} means successful exit.\cr
#'   \code{sys.time} \tab The time elapsed as an object of class \dQuote{proc_time}.
#' }
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
#' \eqn{\alpha_{t} \in R^{m}}{alpha(t) in R^m} \tab
#' \eqn{d_{t} \in R^m}{d(t) in R^m} \tab
#' \eqn{\eta_{t} \in R^m}{eta(t) in R^m} \cr
#' \eqn{T_{t} \in R^{m \times m}}{T(t) in R^{m x m}} \tab
#' \eqn{H_{t} \in R^{m \times m}}{H(t) in R^{m x m}} \tab \cr
#' \eqn{y_{t} \in R^d}{y(t) in R^d} \tab
#' \eqn{c_t \in R^d}{c(t) in R^d} \tab
#' \eqn{\epsilon_{t} \in R^d}{epsilon(t) in R^d} \cr
#' \eqn{Z_{t} \in R^{d \times m}}{Z(t) in R^{d x m}} \tab
#' \eqn{G_{t} \in R^{d \times d}}{G(t) in R^{d x d}} \tab 
#' }
#'
#' Note that \code{fkf} takes as input \code{HHt} and \code{GGt} which
#' corresponds to \eqn{H_t H_t^\prime}{H(t)H(t)'} and \eqn{G_t G_t^\prime}{G(t)G(t)'}. 
#'
#' % <------------------------------------->
#' \strong{Iteration:}
#'
#' The filter iterations are implemented using the expected values
#' \deqn{a_{t} = E[\alpha_t | y_1,\ldots,y_{t-1}]}{a(t) = E[alpha(t) | y(1),...,y(t-1)]}
#' \deqn{a_{t|t} = E[\alpha_t | y_1,\ldots,y_{t}]}{a(t|t) = E[alpha(t) | y(1),...,y(t)]}
#'
#' and variances
#' \deqn{P_{t} = Var[\alpha_t | y_1,\ldots,y_{t-1}]}{P(t) = Var[alpha(t) | y(1),...,y(t-1)]}
#' \deqn{P_{t|t} = Var[\alpha_t | y_1,\ldots,y_{t}]}{P(t|t) = Var[alpha(t) | y(1),...,y(t)]}
#' 
#' of the state \eqn{\alpha_{t}}{alpha(t)} in the following way
#' (for the case of no NA's):
#'
#' Initialisation: Set \eqn{t=1}{t=1} with \eqn{a_{t} = a0}{a(t)=a0} and \eqn{P_{t}=P0}{P(t)=P0}
#'
#' Updating equations:
#' \deqn{v_t = y_t - c_t - Z_t a_t}{v(t) = y(t) - c(t) - Z(t) a(t)}
#' \deqn{F_t = Z_t P_t Z_t^{\prime} + G_t G_t^\prime}{F(t)=Z(t)P(t)Z(t)' + G(t)G(t)'}
#' \deqn{K_t = P_t Z_t^{\prime} F_{t}^{-1}}{K(t) = P(t) Z(t)' F(t)^{-1}}
#' \deqn{a_{t|t} = a_t + K_t v_t}{a(t|t) = a(t) + K(t)v(t)}
#' \deqn{P_{t|t} = P_t - P_t Z_t^\prime K_t^\prime}{P(t|t) = P(t) - P(t) Z(t)' K(t)'}
#'
#' Prediction equations:
#' \deqn{a_{t+1} = d_{t} + T_{t} a_{t|t}}{a(t+1) = d(t) + T(t) a(t|t)}
#' \deqn{P_{t+1} = T_{t} P_{t|t} T_{t}^{\prime} + H_t H_t^\prime}{P(t+1) = T(t)P(t)T(t)' + H(t)H(t)'}
#' 
#' Next iteration: Set \eqn{t=t+1}{t=t+1} and goto \dQuote{Updating equations}.
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

#' 
#' The first element of both \code{at} and \code{Pt} is filled with the
#' function arguments \code{a0} and \code{P0}, and the last, i.e. the (n +
#' 1)-th, element of \code{at} and \code{Pt} contains the predictions for the next time step.
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
#' ## Example: Local level model for the Nile's annual flow.
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
fkf <- function(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt) {

       ## Check all inputs are present
    if(any(c(missing(a0), missing(P0), missing(dt), missing(ct), missing(Tt),
             missing(Zt), missing(HHt), missing(GGt), missing(yt)))){
        stop("None of the input arguments 'a0', 'P0', 'dt', 'ct', 'Tt', 'Zt',",
             "'HHt', 'GGt', and 'yt' must be missing.")
    }

    
    ## Originally the checking of dimensions and type cpuld be disables with 'check.input'
    ## From v 0.17 this options is dropped since it can cause multiple errors as detected by ASAN
    ## From the earlier code:
    ##
    ## 'check.input' should always be 'TRUE' unless the performance
    ## becomes crucial and correctness of the arguments concerning
    ## dimensions, class and storage.mode is ensured.
    ## if( check.input ){
    ##
    ## Now we just issue a warnign about depreciation
    ##if( check.input==FALSE ){
    ##    warning("All inputs are now checked. This flag will be depreciated in a later verion of FKF")
    ##}

 
    ## Check the storage mode: Must be 'double' for all arguments
    stor_mode <- c("a0" = storage.mode(a0), "P0" = storage.mode(P0), "dt" = storage.mode(dt),
                   "ct" = storage.mode(ct), "Tt" = storage.mode(Tt), "Zt" = storage.mode(Zt),
                   "HHt" = storage.mode(HHt), "GGt" = storage.mode(GGt), "yt" = storage.mode(yt))
    
    if(any(stor_mode != "double")){
        stop("Storage mode of variable(s) '",
             paste(names(stor_mode)[stor_mode != "double"],
                   collapse = "', '"),
             "' is not 'double'!\n", sep = "")
    }
    
    ## Check classes of arguments
    error.string <- ""
    if(!is.vector(a0)){
        error.string <- paste(error.string,
                              "'a0' must be of class 'vector'!\n", sep = "")
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
    
    
    time.0 <- proc.time()
    
    
    ans <-.Call("FKF", a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt, PACKAGE = "FKF")
    
    ans$sys.time <- proc.time() - time.0
    ans$yt <- yt
    ans$Zt <- Zt
    ans$Tt <- Tt
    return(ans)
}

