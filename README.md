
<!-- README.md is generated from README.Rmd. Please edit that file -->
MittagLeffleR
=============

The MittagLeffleR `R` package calculates probabilities, quantiles and random variables from both types Mittag-Leffler distributions.

The first type Mittag-Leffler distribution is a heavy-tailed distribution, and occurs mainly as a waiting time distribution in problems with "fractional" time scales, e.g. times between earthquakes.

The second type Mittag-Leffler distribution is light-tailed and, in a sense, inverse to the family of sum-stable distributions. It is used for time-changes of stochastic processes, to create "time-fractional" evolutions, e.g. anomalous diffusion processes.

Installation
------------

You can install MittagLeffleR from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("strakaps/MittagLeffler")
```

Examples
--------

### Fitting a Mittag-Leffler distribution to a random dataset

Fit a Mittag-Leffler distribution (first type) to a collection of Mittag-Leffler random variables:

``` r
library(MittagLeffleR)
y = rml(n = 100, tail = 0.9, scale = 2)
mlml <- function(X) {
  log_l <- function(theta) {
    #transform parameters so can do optimization unconstrained
    theta[1] <- 1/(1+exp(-theta[1]))
    theta[2] <- exp(theta[2])
    -sum(log(dml(X,theta[1],theta[2])))
  }
  ml_theta <- stats::optim(c(0.5,0.5), fn=log_l)$par
  #transform back
  ml_a <- 1/(1 + exp(-ml_theta[1]))
  ml_d <- exp(ml_theta[2])
  return(list("tail" = ml_a, "scale" = ml_d))
}
mlml(y)
#> $tail
#> [1] 0.8924838
#> 
#> $scale
#> [1] 2.370469
```

### Calculate the probability density of an anomalous diffusion process

Brownian motion at time *t* has a normal probability density *n*(*x*, *t*<sup>1/2</sup>). A fractional diffusion at time *t* has the time-changed probability density

*p*(*x*, *t*)=∫*n*(*x*, *u*)*h*(*u*, *t*)*d**u*

where *h*(*u*, *t*) is a second type Mittag-Leffler probability density with scale *t*<sup>*α*</sup>:

``` r
tail <- 0.65
dx <- 0.01
x <- seq(-3,3,dx)
umax <- qml(p = 0.99, tail = tail, scale = 1, second.type = TRUE)
u <- seq(0.01,umax,dx)
h <- dml(x = u, tail = tail, scale = 1, second.type = TRUE)
N <- outer(x,u,function(x,u){dnorm(x = x, sd = sqrt(u))})
p <- N %*% h * dx
plot(x,p, type='l', main = "Density of fractional diffusion at t=1")
```

![](README-unnamed-chunk-3-1.png)
