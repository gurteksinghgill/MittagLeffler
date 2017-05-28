#' Mittag-Leffler Distribution
#'
#' These functions provide probability density, cumulative distribution
#' function, quantile function and random variate generation for the 
#' two types of Mittag-Leffler distribution.
#'
#' @param t vector of quantiles
#' @param p vector of probabilities
#' @param a alpha parameter of Mittag-Leffler distribution
#' @param d delta (scale) parameter of Mittag-Leffler distribution
#' @param n Number of random values
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE, probabilities are P(X ≤ x) 
#'     otherwise, P(X > x)
#' @details
#' The generalized Mittag-Leffer function is defined by the power series 
#'     \deqn{E_(\alpha,\beta) (z) = \sum_(k=0)^(\inf) (z^k)/\Gamma(\alpha 
#'     k + \beta) }
#' for complex \eqn{z} and complex \eqn{\alpha, \beta} with 
#' \eqn{Real(\alpha) > 0}.
#'
#' The Mittag-Leffler CDF is given by 
#' \deqn{F(t; \alpha, \delta) = 1 - E_(\alpha,1) (-(t/\delta)^\alpha)}
#' for \eqn{t \ge 0}, \eqn{0 < \alpha \le 1} 
#' and scale parameter \eqn{\delta > 0}.
#'
#' The Mittag-Leffler PDF is given by 
#' \deqn{f(t; \alpha, \delta) = (t^(\alpha - 1))/(\delta^\alpha) 
#' E_(\alpha,\alpha) [-(t/\delta)^\alpha]}
#' for \eqn{t \ge 0}, \eqn{0 < \alpha \le 1} 
#' and scale parameter \eqn{\delta > 0}.
#' @return \code{dml} returns the density, 
#'         \code{pml} returns the distribution function, 
#'         \code{qml} returns the quantile function, and 
#'         \code{rml} generates random values.
#' @source Uses Garrapa's code...
#' @references Garrapa's paper and others...
#' @seealso Anything here?
#' @name MittagLeffler
NULL


#' @rdname MittagLeffler
#' @examples
#' dml(1, 0.8)
#' @export
dml <- function(t,a,d=1,log=FALSE) {
  ml <- (t^(a-1)/(d^a))*mlf(-(t/d)^a, a, a, 1)
  if(log==FALSE) {
    return(ml)
  }
  else {
    return(log(ml))
  }
}



#' @rdname MittagLeffler
#' @examples
#' pml(2, 0.7, 1.5)
#' @export
pml <- function(t,a,d=1,lower.tail=TRUE, log.p=FALSE) {
  ml <- 1 - mlf(-(t/d)^a,a,1,1)
  if(lower.tail==FALSE) {ml <- 1 - ml}
  if(log.p==FALSE) {
    return(ml)
  }
  else {
    return(log(ml))
  }
}



#' @rdname MittagLeffler
#' @examples
#' qml(0.25, 0.9)
#' @export
qml <- function(p,a,d=1,lower.tail=TRUE, log.p=FALSE) {
  if(log.p==TRUE) {
    p <- exp(p)
  }
  if(lower.tail==FALSE) {p <- 1 - p}
  x <- numeric(length(p))
  for (i in 1:length(p)) {
    qml_p <- function(t) {pml(t,a,d) - p[i]}
    x[i] <- stats::uniroot(qml_p, interval = c(10^-14,100), extendInt="upX")$root
  }
  return(x)
}



#' @rdname MittagLeffler
#' @examples
#' rml(1000, 0.7, 1)
#'
#' ##Approximating Mittag-Leffler distribution parameters
#' ##alpha and delta for observations X by Maximum Likelihood
#'
#' mlml <- function(X) {
#'   log_l <- function(theta) {
#'     #transform parameters so can do optimization unconstrained
#'     theta[1] <- 1/(1+exp(-theta[1]))
#'     theta[2] <- exp(theta[2])
#'     -sum(log(dml(X,theta[1],theta[2])))
#'   }
#'   ml_theta <- stats::optim(c(0.5,0.5), fn=log_l)$par
#'   ml_a <- 1/(1 + exp(-ml_theta[1]))
#'   ml_d <- exp(ml_theta[2])
#'   print(paste("alpha =", ml_a, "delta =", ml_d))
#' }
#' @export
rml <- function(n,a,d=1){
  x <- numeric(n)
  u <- stats::runif(n)
  for (i in 1:n) {
    pml_u <- function(t) {pml(t,a,d) - u[i]}
    x[i] <- stats::uniroot(pml_u, interval = c(10^-14,100), extendInt="upX")$root
  }
  return(x)
}

