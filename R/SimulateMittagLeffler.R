#' The Mittag-Leffler Distribution
#'
#' These functions provide density, distribution function, quantile function and random generation for the Mittag-Leffler distribution
#' with parameters alpha and delta.
#'
#' @param t vector of quantiles
#' @param p vector of probabilities
#' @param a alpha parameter of Mittag-Leffler distribution
#' @param d delta (scale) parameter of Mittag-Leffler distribution
#' @param n Number of random values
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE, probabilities are P(X ≤ x) otherwise, P(X > x)
#' @details
#' The generalized Mittag-Leffer function is defined by the power series \deqn{E_(\alpha,\beta) (z) = \sum_(k=0)^(\inf) (z^k)/\Gamma(\alpha k + \beta) }
#' for complex \eqn{z} and complex \eqn{\alpha, \beta} with \eqn{Real(\alpha) > 0}.
#'
#' The Mittag-Leffler CDF is given by \deqn{F(t; \alpha, \delta) = 1 - E_(\alpha,1) (-(t/\delta)^\alpha)}
#' for \eqn{t \ge 0}, \eqn{0 < \alpha \le 1} and scale parameter \eqn{\delta > 0}.
#'
#' The Mittag-Leffler PDF is given by \deqn{f(t; \alpha, \delta) = (t^(\alpha - 1))/(\delta^\alpha) E_(\alpha,\alpha) [-(t/\delta)^\alpha]}
#' for \eqn{t \ge 0}, \eqn{0 < \alpha \le 1} and scale parameter \eqn{\delta > 0}.
#' @return \code{dml} returns the density, \code{pml} returns the distribution function, \code{qml} returns the quantile function, and \code{rml} generates random values.
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












mlf <- function(z,a,b=1,g=1) {
  #check number of inputs and handle errors
  if (nargs() < 2) {
    #handle error
  }

  #check inputs values and handle errors

  #target precision
  log_epsilon <- log(10^(-15))

  #LT inversion for each element of z
  E <- numeric(length(z))
  for (i in 1:length(z)) {
    if (abs(z[i]) < 10^(-15) ) {
      E[i] <- 1/gamma(b)
    }
    else {
      E[i] <- LTinversion(1,z[i],a,b,g,log_epsilon)
    }
  }

  return(E)
}

LTinversion <- function(t, lambda, a, b, g, log_epsilon) {
  theta <- Arg(lambda)
  kmin <- ceiling(-a/2 - theta/2/pi)
  kmax <- floor(a/2 - theta/2/pi)
  if (kmin > kmax) {
    k_vett <- c()
  }
  else {
  k_vett <- kmin:kmax
  }
  s_star <- abs(lambda)^(1/a) * exp(1i*(theta+2*k_vett*pi)/a)

  phi_s_star <- (Re(s_star)+abs(s_star))/2

  phi_s_star <- sort(phi_s_star, index.return=TRUE)$x
  index_s_star <- sort(phi_s_star, index.return=TRUE)$ix
  s_star <- s_star[index_s_star]

  index_save <- phi_s_star > 10^(-15)
  s_star <- s_star[index_save]
  phi_s_star <- phi_s_star[index_save]

  s_star <- c(0, s_star)
  phi_s_star <- c(0, phi_s_star)
  J1 <- length(s_star)
  J <- J1 - 1

  p <- c( max(0,-2*(a*g-b+1)) , rep(1,J)*g )
  q <- c( rep(1,J)*g , +Inf)
  phi_s_star <- c(phi_s_star, +Inf)

  admissible_regions <- which( (phi_s_star[1:length(phi_s_star)-1] < (log_epsilon - log(.Machine$double.eps))/t) & (phi_s_star[1:length(phi_s_star)-1] < phi_s_star[2:length(phi_s_star)]))

  JJ1 <- admissible_regions[length(admissible_regions)]
  mu_vett <- rep(1,JJ1)*Inf
  N_vett <- rep(1,JJ1)*Inf
  h_vett <- rep(1,JJ1)*Inf

  find_region <- 0
  while (!find_region) {
    for (j1 in admissible_regions) {
        if (j1 < J1) {
            temp2 <- OptimalParam_RB(t,phi_s_star[j1],phi_s_star[j1+1],p[j1],q[j1],log_epsilon)
            muj <- temp2[1]
            hj <- temp2[2]
            Nj <- temp2[3]
        }
        else {
            temp3 <- OptimalParam_RU(t,phi_s_star[j1],p[j1],log_epsilon)
            muj <- temp3[1]
            hj <- temp3[2]
            Nj <- temp3[3]
        }
        mu_vett[j1] <- muj
        h_vett[j1] <- hj
        N_vett[j1] <- Nj
    }
    if (min(N_vett) > 200) {
        log_epsilon <- log_epsilon +log(10)
    }
    else {
        find_region <- 1
    }
  }

  N <- min(N_vett)
  iN <- which.min(N_vett)
  mu <- mu_vett[iN]
  h <- h_vett[iN]

  k <- -N:N
  u <- h*k
  z <- mu*(1i*u+1)^2
  zd <- -2*mu*u + 2*mu*1i
  zexp <- exp(z*t)
  F <- z^(a*g-b)/(z^a - lambda)^g*zd
  S <- zexp*F
  Integral <- h*sum(S)/2/pi/1i

  if ((iN+1) > length(s_star)) {
    ss_star <- c()
  }
  else {
    ss_star <- s_star[(iN+1):length(s_star)]
  }
  Residues <- sum(1/a*(ss_star)^(1-b)*exp(t*ss_star))

  E <- Integral + Residues
  if (!is.complex(lambda)) {
    E <- Re(E)
  }
  return(E)
}

OptimalParam_RB <- function(t, phi_s_star_j, phi_s_star_j1, pj, qj, log_epsilon) {
  log_eps <- -36.043653389117154 #log(eps)
  fac <- 1.01
  conservative_error_analysis <- 0

  f_max <- exp(log_epsilon - log_eps)

  sq_phi_star_j <- sqrt(phi_s_star_j)
  threshold <- 2*sqrt((log_epsilon - log_eps)/t) ;
  sq_phi_star_j1 <- min(sqrt(phi_s_star_j1), threshold - sq_phi_star_j) ;

  if ((pj < 10^(-14)) && (qj < 10^(-14))) {
    sq_phibar_star_j <- sq_phi_star_j
    sq_phibar_star_j1 <- sq_phi_star_j1
    adm_region <- 1
  }

  if ((pj < 10^(-14)) && (qj >= 10^(-14))) {
    sq_phibar_star_j <- sq_phi_star_j
    if (sq_phi_star_j > 0) {
        f_min <- fac*(sq_phi_star_j/(sq_phi_star_j1-sq_phi_star_j))^qj
    }
    else {
        f_min <- fac
    }
    if (f_min < f_max) {
        f_bar <- f_min + f_min/f_max*(f_max-f_min)
        fq <- f_bar^(-1/qj)
        sq_phibar_star_j1 <- (2*sq_phi_star_j1 - fq*sq_phi_star_j)/(2+fq)
        adm_region <- 1
    }
    else {
        adm_region <- 0
    }
  }

  if ((pj >= 10^(-14)) && (qj < 10^(-14))) {
    sq_phibar_star_j1 <- sq_phi_star_j1
    f_min = fac*(sq_phi_star_j1/(sq_phi_star_j1 - sq_phi_star_j))^pj
    if (f_min < f_max) {
        f_bar <- f_min + f_min/f_max*(f_max - f_min)
        fp <- f_bar^(-1/pj)
        sq_phibar_star_j <- (2*sq_phi_star_j+fp*sq_phi_star_j1)/(2-fp)
        adm_region <- 1
    }
    else {
        adm_region <- 0
    }
  }

  if ((pj >= 10^(-14)) && (qj >= 10^(-14))) {
    f_min = fac*(sq_phi_star_j+sq_phi_star_j1)/(sq_phi_star_j1-sq_phi_star_j)^max(pj,qj)
    if (f_min < f_max) {
        f_min <- max(f_min,1.5)
        f_bar <- f_min + f_min/f_max*(f_max-f_min)
        fp <- f_bar^(-1/pj)
        fq <- f_bar^(-1/qj)
        if (!conservative_error_analysis) {
            w <- -phi_s_star_j1*t/log_epsilon
        }
        else {
            w <- -2*phi_s_star_j1*t/(log_epsilon-phi_s_star_j1*t)
        }
        den <- 2+w - (1+w)*fp + fq
        sq_phibar_star_j <- ((2+w+fq)*sq_phi_star_j + fp*sq_phi_star_j1)/den
        sq_phibar_star_j1 <- (-(1+w)*fq*sq_phi_star_j + (2+w-(1+w)*fp)*sq_phi_star_j1)/den
        adm_region <- 1
    }
    else {
        adm_region <- 0
    }
  }

  if (adm_region) {
    log_epsilon <- log_epsilon  - log(f_bar)
    if (!conservative_error_analysis) {
        w <- -sq_phibar_star_j1^2*t/log_epsilon
    }
    else {
        w = -2*sq_phibar_star_j1^2*t/(log_epsilon-sq_phibar_star_j1^2*t)
    }
    muj <- (((1+w)*sq_phibar_star_j + sq_phibar_star_j1)/(2+w))^2
    hj <- -2*pi/log_epsilon*(sq_phibar_star_j1-sq_phibar_star_j)/((1+w)*sq_phibar_star_j + sq_phibar_star_j1)
    Nj <- ceiling(sqrt(1-log_epsilon/t/muj)/hj)
  }
  else {
    muj <- 0
    hj <- 0
    Nj <- +Inf
  }

  out <- c(muj, hj, Nj)
  return(out)

}

OptimalParam_RU <- function(t, phi_s_star_j, pj, log_epsilon) {
  sq_phi_s_star_j <- sqrt(phi_s_star_j)
  if (phi_s_star_j > 0) {
    phibar_star_j <- phi_s_star_j*1.01
  }
  else {
    phibar_star_j <- 0.01
  }
  sq_phibar_star_j <- sqrt(phibar_star_j)

  f_min <- 1
  f_max <- 10
  f_tar <- 5

  stopp <- 0
  while (!stopp) {
    phi_t <- phibar_star_j*t
    log_eps_phi_t <- log_epsilon/phi_t
    Nj <- ceiling(phi_t/pi*(1 - 3*log_eps_phi_t/2 + sqrt(1-2*log_eps_phi_t)))
    A <- pi*Nj/phi_t
    sq_muj <- sq_phibar_star_j*abs(4-A)/abs(7-sqrt(1+12*A))
    fbar <- ((sq_phibar_star_j-sq_phi_s_star_j)/sq_muj)^(-pj)
    stopp <- (pj < 10^(-14)) || (f_min < fbar && fbar < f_max)
    if (!stopp) {
        sq_phibar_star_j <- f_tar^(-1/pj)*sq_muj + sq_phi_s_star_j
        phibar_star_j <- sq_phibar_star_j^2
    }
  }
  muj <- sq_muj^2
  hj <- (-3*A - 2 + 2*sqrt(1+12*A))/(4-A)/Nj

  log_eps <- log(.Machine$double.eps)
  threshold <- (log_epsilon - log_eps)/t
  if (muj > threshold) {
    if (abs(pj) < 10^(-14)) {
      Q <- 0
    }
    else {
      Q <- f_tar^(-1/pj)*sqrt(muj)
    }
    phibar_star_j <- (Q + sqrt(phi_s_star_j))^2
    if (phibar_star_j < threshold) {
        w <- sqrt(log_eps/(log_eps-log_epsilon))
        u <- sqrt(-phibar_star_j*t/log_eps)
        muj <- threshold
        Nj <- ceiling(w*log_epsilon/2/pi/(u*w-1))
        hj <- sqrt(log_eps/(log_eps - log_epsilon))/Nj
    }
    else {
        Nj <- +Inf
        hj = 0
    }
  }

  out <- c(muj, hj, Nj)
  return(out)
}


