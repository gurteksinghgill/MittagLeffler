#' Mittag-Leffler Distribution
#'
#' Probability density, cumulative distribution
#' function, quantile function and random variate generation for the
#' two types of Mittag-Leffler distribution.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken
#'        to be the number required.
#' @param alpha tail parameter.
#' @param tau scale parameter.
#' @param second.type logical; if FALSE (default), 
#'        first type of Mittag-Leffler distribution is assumed.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities are P[X ≤ x]
#'        otherwise, P[X > x]

#' @details
#' The generalized (two-parameter) Mittag-Leffer function is defined by the
#' power series
#'     \deqn{E_{\alpha,\beta} (z) = \sum_{k=0}^\infty  z^k / \Gamma(\alpha
#'     k + \beta) }
#' for complex \eqn{z} and complex \eqn{\alpha, \beta} with
#' \eqn{Real(\alpha) > 0}.
#'
#' The **first type** of Mittag-Leffler distribution assumes the Mittag-Leffler
#' function as its tail function, so that the CDF is given by
#' \deqn{F(t; \alpha, \tau) = 1 - E_{\alpha,1} (-(t/\tau)^\alpha)}
#' for \eqn{t \ge 0}, \eqn{0 < \alpha \le 1}
#' and scale parameter \eqn{\tau > 0}.
#' Its PDF is given by
#' \deqn{f(t; \alpha, \tau) = t^{\alpha - 1} 
#' E_{\alpha,\alpha} [-(t/\tau)^\alpha] / \tau^\alpha.}
#' As \eqn{\alpha} approaches 1 from below, the Mittag-Leffler converges
#' (weakly) to the expnoential
#' distribution. For \eqn{0 < \alpha < 1}, it is (very) heavy-tailed, i.e.
#' has infinite mean.
#'
#' The **second type** of Mittag-Leffler distribution is defined via the
#' Laplace transform of its density f:
#' \deqn{\int_0^\infty \exp(-st) f(t; \alpha, 1) dt = E_{\alpha,1}(-s)}
#' It is light-tailed, i.e. all its moments are finite.
#'
#' @return \code{dml} returns the density,
#'         \code{pml} returns the distribution function,
#'         \code{qml} returns the quantile function, and
#'         \code{rml} generates random variables.

#' @references
#' Haubold, H. J., Mathai, A. M., & Saxena, R. K. (2011). Mittag-Leffler
#' Functions and Their Applications. Journal of Applied Mathematics, 2011, 
#' 1–51. \url{http://doi.org/10.1155/2011/298628}
#' 
#' Garrappa, R. (2015). Numerical Evaluation of Two and Three Parameter
#' Mittag-Leffler Functions. SIAM Journal on Numerical Analysis, 53(3),
#'  1350–1369. \url{http://doi.org/10.1137/140971191}
#'
#' Mittag-Leffler distribution. (2017, May 3).
#' In Wikipedia, The Free Encyclopedia.
#' \url{https://en.wikipedia.org/w/index.php?title=Mittag-Leffler_distribution&oldid=778429885}
#'
#' The Mittag-Leffler function. MathWorks File Exchange.
#' \url{https://au.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function}

#' @name MittagLeffler
NULL


#' @rdname MittagLeffler
#' @examples
#' dml(1, 0.8)
#' @export
dml <- function(x,alpha,scale=1,log=FALSE, second.type=FALSE){
  if (length(alpha) > 1){
    stop("length(alpha) must be 1.")
  }
  if (second.type==FALSE) {
    return(dml1(x,alpha,scale,log))
  } else {
    return(dml2(x,alpha,scale,log))
  }
}

# first type
dml1 <- function(t,alpha,scale=1,log=FALSE) {
  ml <- (t^(alpha-1)/(scale^alpha))*mlf(-(t/scale)^alpha, alpha, alpha, 1)
  if (log) {
    ml <- log(ml)
  }
  return(ml)
}

# second type
dml2 <- function(u,alpha,scale=1,log=FALSE) {
  # find the distribution of E(scale^(1/alpha)), where E() is the inverse stable
  # subordinator; Meerschaert & Straka, Eq.(9)
  t = scale^(1/alpha)
  h=(t/alpha)*u^(-1-1/alpha)*stabledist::dstable(t*u^(-1/alpha), alpha=alpha, beta=1.0,
                             gamma=1.0, delta=0.0, pm=1)
  if (log) {
    h <- log(h)
  }
  return(h)
}


#' @rdname MittagLeffler
#' @examples
#' pml(2, 0.7, 1.5)
#' @export
pml <- function(q, alpha, scale=1, second.type=FALSE, lower.tail=TRUE, 
                log.p=FALSE) {
  # rescale
  q <- q/scale
  if (!second.type){
    p <- pml1(q,alpha)
  } else {
    p <- pml2(q,alpha)
  }
  if (!lower.tail) {
    p <- 1-p
  }
  if (log.p) {
    p <- log(p)
  }
  return(p)
}

# type 1 with unit scale
pml1 <- function(q,alpha) {
  p <- 1 - mlf(-(q)^alpha,alpha,1,1)
}

# type 2 with unit scale
pml2 <- function(q,alpha) {
  p <- stabledist::pstable(q^(-1/alpha), alpha=alpha, beta=1, gamma=1, delta=0,
                           pm=1, lower.tail = FALSE)
}

#' @rdname MittagLeffler
#' @examples
#' qml(p = c(0.25, 0.5, 0.75), alpha = 0.6, scale = 100)
#' @export

qml <- function(p, alpha, scale=1, second.type=FALSE, lower.tail=TRUE,
                log.p=FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1-p
  }
  if (!second.type) {
    q <- qml1(p, alpha, scale)
  } else {
    q <- scale * qml2(p, alpha)
  }
  return(q)
}

qml1 <- function(p,alpha,scale=1) {
  x <- numeric(length(p))
  for (i in 1:length(p)) {
    qml_p <- function(t) {pml(t,alpha,scale) - p[i]}
    x[i] <- stats::uniroot(qml_p, interval = c(10^-14,100),
                           extendInt="upX")$root
  }
  return(x)
}

# type 2 with unit scale
qml2 <- function(p, alpha){
  q <- stabledist::qstable(p, alpha=alpha, beta=1, gamma=1, delta=0, pm=1,
                           lower.tail=FALSE)^(-alpha)
}

#' @rdname MittagLeffler
#' @examples
#' rml(1000, 0.7, 1)
#'
#' ## Fitting Mittag-Leffler distribution to observations X by Maximum
#' ## Likelihood
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
#'   print(paste("alpha =", ml_a, "scale =", ml_d))
#' }
#' mlml(rml(n = 100, alpha = 0.9, scale = 2))

#' @export
rml <- function(n,alpha,scale=1, second.type=FALSE){
  if (length(n) > 1){
    n <- length(n)
  }
  if (!second.type){
    x <- scale * rml1(n,alpha)
  } else {
    x <- scale * rml2(n,alpha)
  }
  return(x)
}

# unit scale; see e.g. Haubold, Mathai & Saxena (2011)
rml1 <- function(n, alpha){
  gamma <- (cos(pi*alpha/2))^(1/alpha)
  y <- stabledist::rstable(n, alpha=alpha, beta=1, gamma=gamma, delta=0, pm=1)
  x <- stats::rexp(n)
  y * x^(1/alpha)
}

rml2 <- function(n, alpha){
  stabledist::rstable(n, alpha=alpha, beta=1, gamma=1, delta=0, pm=1)^(-alpha)
}
