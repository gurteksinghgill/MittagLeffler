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
#' @param second.type logical; if FALSE, first type of Mittag-Leffler
#'     distribution is assumed
#' @details
#' The generalized (two-parameter) Mittag-Leffer function is defined by the
#' power series 
#'     \deqn{E_(\alpha,\beta) (z) = \sum_(k=0)^(\inf) (z^k)/\Gamma(\alpha 
#'     k + \beta) }
#' for complex \eqn{z} and complex \eqn{\alpha, \beta} with 
#' \eqn{Real(\alpha) > 0}.
#'
#' The **first type** of Mittag-Leffler distribution assumes the Mittag-Leffler
#' function as its tail function, so that the CDF is given by 
#' \deqn{F(t; \alpha, \delta) = 1 - E_(\alpha,1) (-(t/\delta)^\alpha)}
#' for \eqn{t \ge 0}, \eqn{0 < \alpha \le 1} 
#' and scale parameter \eqn{\delta > 0}.
#' Its PDF is given by 
#' \deqn{f(t; \alpha, \delta) = (t^(\alpha - 1))/(\delta^\alpha) 
#' E_(\alpha,\alpha) [-(t/\delta)^\alpha].}
#' For \eqn{\alpha = 1}, the Mittag-Leffler is identical to the expnoential
#' distribution. For \eqn{0 < \alpha < 1}, it is (very) heavy-tailed, i.e.
#' has infinite mean. 
#' 
#' The **second type** of Mittag-Leffler distribution is defined via its
#' Laplace transform: 
#' \deqn{\int_0^{\infty} e^{-st} f(t; \alpha, 1) dt = E_{\alpha,1}(-s)}
#' It is light-tailed, i.e. all its moments are finite. 
#' 
#' @return \code{dml} returns the density, 
#'         \code{pml} returns the distribution function, 
#'         \code{qml} returns the quantile function, and 
#'         \code{rml} generates random values.
#' @references
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
dml <- function(q,a,d=1,p.log=FALSE, second.type=FALSE){
  if (second.type==FALSE) {
    return(dml1(q,a,d,p.log))
  } else {
    return(dml2(q,a,d,p.log))
  }
}

# first type
dml1 <- function(t,a,d=1,p.log=FALSE) {
  ml <- (t^(a-1)/(d^a))*mlf(-(t/d)^a, a, a, 1)
  if (p.log) {
    ml <- log(ml)
  }
  return(ml)
}

# second type
dml2 <- function(u,a,d=1,p.log=FALSE) {
  # find the distribution of E(d^(1/a)), where E() is the inverse stable
  # subordinator; Meerschaert & Straka, Eq.(9)
  t = d^(1/a)
  h=(t/a)*u^(-1-1/a)*stabledist::dstable(t*u^(-1/a), alpha=a, beta=1.0, 
                             gamma=1.0, delta=0.0, pm=1)
  if (p.log) {
    h <- log(h)
  }
  return(h)
}


#' @rdname MittagLeffler
#' @examples
#' pml(2, 0.7, 1.5)
#' @export
pml <- function(q,a,d=1, second.type=FALSE, lower.tail=TRUE, log.p=FALSE) {
  # rescale
  q <- q/d
  if (!second.type){
    p <- pml1(q,a)
  } else {
    p <- pml2(q,a)
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
pml1 <- function(q,a) {
  p <- 1 - mlf(-(q)^a,a,1,1)
}

# type 2 with unit scale
pml2 <- function(q,a) {
  p <- stabledist::pstable(q^(-1/a), alpha=a, beta=1, gamma=1, delta=0, 
                           pm=1, lower.tail = FALSE)
}

#' @rdname MittagLeffler
#' @examples
#' qml(0.25, 0.9)
#' @export

qml <- function(p, a, d=1, second.type=FALSE, lower.tail=TRUE, 
                log.p=FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1-p
  }
  if (!second.type) {
    q <- qml1(p, a, d)
  } else {
    q <- d * qml2(p, a)
  }
  return(q)
}

qml1 <- function(p,a,d=1) {
  x <- numeric(length(p))
  for (i in 1:length(p)) {
    qml_p <- function(t) {pml(t,a,d) - p[i]}
    x[i] <- stats::uniroot(qml_p, interval = c(10^-14,100), 
                           extendInt="upX")$root
  }
  return(x)
}

# type 2 with unit scale
qml2 <- function(p, a){
  q <- stabledist::qstable(p, alpha=a, beta=1, gamma=1, delta=0, pm=1,
                           lower.tail=FALSE)^(-a)
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
rml <- function(n,a,d=1, second.type=FALSE){
  if (!second.type){
    x <- d * rml1(n,a)
  } else {
    x <- d * rml2(n,a)
  }
  return(x)
}

# unit scale; see e.g. Haubold, Mathai & Saxena (2011)
rml1 <- function(n, a){
  scale <- (cos(pi*a/2))^(1/a)
  y <- stabledist::rstable(n, alpha=a, beta=1, gamma=scale, delta=0, pm=1)
  x <- stats::rexp(n)
  y * x^(1/a)
}

rml2 <- function(n, a){
  stabledist::rstable(n, alpha=a, beta=1, gamma=1, delta=0, pm=1)^(-a)
}