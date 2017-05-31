library(MittagLeffler)
context("Does dml() integrate to pml()?")

alphavec <- seq(0.4, 0.9, 0.1)
scalevec <- 10^seq(-3,3,2)
pvec <- seq(0.1, 0.9, 0.1)
tol <- 0.005

test_that("dml() integrates to pml() for Type 1", {
  for (alpha in alphavec){
    for (scale in scalevec){
      ml_dens <- function(x) {
        dml(x = x, alpha = alpha, scale = scale, second.type = FALSE)
      }
      qvec <- qml(pvec, alpha, scale, second.type = FALSE)
      for (i in 1:(length(qvec)-1)){
                p2 <- integrate(ml_dens, qvec[i], qvec[i+1])$value
        expect_equal(object = p2, expected = pvec[i+1]-pvec[i], tol=tol)
      }
    }
  }
})

test_that("dml() integrates to pml() for Type 2", {
  for (alpha in alphavec){
    for (scale in scalevec){
      qvec <- qml(pvec, alpha, scale, second.type = TRUE)
      ml_dens <- function(x) {
        dml(x = x, alpha = alpha, scale = scale, second.type = TRUE)
      }
      for (i in 1:(length(qvec)-1)){
        p2 <- integrate(ml_dens, qvec[i], qvec[i+1])$value
        expect_equal(object = p2, expected = pvec[i+1]-pvec[i], tol=tol)
      }
    }
  }
})

