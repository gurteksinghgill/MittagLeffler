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
#' @param tail tail parameter.
#' @param scale scale parameter.
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
#' \deqn{F(q; \alpha, \tau) = 1 - E_{\alpha,1} (-(q/\tau)^\alpha)}
#' for \eqn{q \ge 0}, tail parameter \eqn{0 < \alpha \le 1},
#' and scale parameter \eqn{\tau > 0}.
#' Its PDF is given by
#' \deqn{f(x; \alpha, \tau) = x^{\alpha - 1} 
#' E_{\alpha,\alpha} [-(x/\tau)^\alpha] / \tau^\alpha.}
#' As \eqn{\alpha} approaches 1 from below, the Mittag-Leffler converges
#' (weakly) to the expnoential
#' distribution. For \eqn{0 < \alpha < 1}, it is (very) heavy-tailed, i.e.
#' has infinite mean.
#'
#' The **second type** of Mittag-Leffler distribution is defined via the
#' Laplace transform of its density f:
#' \deqn{\int_0^\infty \exp(-sx) f(x; \alpha, 1) dx = E_{\alpha,1}(-s)}
#' It is light-tailed, i.e. all its moments are finite.
#' At scale \eqn{\tau}, its density is 
#' \deqn{f(x; \alpha, \tau) = f(x/\tau; \alpha, 1) / \tau.}
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
#' dml(1, 0.6, second.type=TRUE)
#' @export
dml <- function(x,tail,scale=1,log=FALSE, second.type=FALSE){
  if (length(tail) > 1){
    stop("length(tail) must be 1.")
  }
  if (second.type==FALSE) {
    return(dml1(x,tail,scale,log))
  } else {
    y <- dml2(x/scale,tail)/scale
    if (!log) return(y) else return(log(y))
  }
}

# first type
dml1 <- function(t,tail,scale=1,log=FALSE) {
  ml <- (t^(tail-1)/(scale^tail))*mlf(-(t/scale)^tail, tail, tail, 1)
  if (log) {
    ml <- log(ml)
  }
  return(ml)
}

# second type, unit scale
dml2 <- function(u,tail) {
  # find the distribution of E(1), where E() is the inverse stable
  # subordinator; Meerschaert & Straka, Eq.(9).
  # The scale parameter in the Samorodnitsky & Taqqu representation which
  # makes the stable distribution have Laplace transform exp(-s^tail):
  gamma <- (cos(pi*tail/2))^(1/tail)
  
  h=(1/tail)*u^(-1-1/tail)*
  #  stable::dstable(x = u^(-1/tail), loc = 0, disp = 1, skew = 1, 
  #                  tail = tail, eps = 1.0e-6)
     stabledist::dstable(x = u^(-1/tail), alpha = tail, beta = 1, 
                         gamma = gamma, delta = 0, pm = 1)
}


#' @rdname MittagLeffler
#' @examples
#' pml(2, 0.7, 1.5)
#' @export
pml <- function(q, tail, scale=1, second.type=FALSE, lower.tail=TRUE, 
                log.p=FALSE) {
  # rescale
  q <- q/scale
  if (!second.type){
    p <- pml1(q,tail)
  } else {
    p <- pml2(q,tail)
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
pml1 <- function(q,tail) {
  p <- 1 - mlf(-(q)^tail,tail,1,1)
}

# type 2 with unit scale
pml2 <- function(q,tail) {
  # the scale parameter in the Samorodnitsky & Taqqu representation which
  # makes the stable distribution have Laplace transform exp(-s^tail):
  gamma <- (cos(pi*tail/2))^(1/tail)
  p <- stabledist::pstable(q^(-1/tail), alpha=tail, beta=1, gamma=gamma, delta=0,
                           pm=1, lower.tail = FALSE)
}

#' @rdname MittagLeffler
#' @examples
#' qml(p = c(0.25, 0.5, 0.75), tail = 0.6, scale = 100)
#' @export

qml <- function(p, tail, scale=1, second.type=FALSE, lower.tail=TRUE,
                log.p=FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1-p
  }
  if (!second.type) {
    q <- qml1(p, tail, scale)
  } else {
    q <- scale * qml2(p, tail)
  }
  return(q)
}

qml1 <- function(p,tail,scale=1) {
  x <- numeric(length(p))
  for (i in 1:length(p)) {
    qml_p <- function(t) {pml(t,tail,scale) - p[i]}
    x[i] <- stats::uniroot(qml_p, interval = c(10^-14,100),
                           extendInt="upX", tol = 1e-14)$root
  }
  return(x)
}

# type 2 with unit scale
qml2 <- function(p, tail){
  # the scale parameter in the Samorodnitsky & Taqqu representation which
  # makes the stable distribution have Laplace transform exp(-s^tail):
  gamma <- (cos(pi*tail/2))^(1/tail)
  q <- stabledist::qstable(p, alpha=tail, beta=1, gamma=gamma, delta=0, pm=1,
                           lower.tail=FALSE)^(-tail)
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
#'   print(paste("tail =", ml_a, "scale =", ml_d))
#' }
#' mlml(rml(n = 100, tail = 0.9, scale = 2))

#' @export
rml <- function(n,tail,scale=1, second.type=FALSE){
  if (length(n) > 1){
    n <- length(n)
  }
  if (!second.type){
    x <- scale * rml1(n,tail)
  } else {
    x <- scale * rml2(n,tail)
  }
  return(x)
}

# unit scale; see e.g. Haubold, Mathai & Saxena (2011)
rml1 <- function(n, tail){
  # the scale parameter in the Samorodnitsky & Taqqu representation which
  # makes the stable distribution have Laplace transform exp(-s^tail):
  gamma <- (cos(pi*tail/2))^(1/tail)
  y <- stabledist::rstable(n, alpha=tail, beta=1, gamma=gamma, delta=0, pm=1)
  x <- stats::rexp(n)
  y * x^(1/tail)
}

rml2 <- function(n, tail){
  # the scale parameter in the Samorodnitsky & Taqqu representation which
  # makes the stable distribution have Laplace transform exp(-s^tail):
  gamma <- (cos(pi*tail/2))^(1/tail)
  stabledist::rstable(n, alpha=tail, beta=1, gamma=gamma, delta=0, pm=1)^(-tail)
}
