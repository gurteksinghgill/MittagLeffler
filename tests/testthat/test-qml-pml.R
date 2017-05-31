library(MittagLeffler)
context("Compatibility of qml() and pml()")

alphavec <- seq(0.4, 0.9, 0.1)
scalevec <- 10^seq(-3,3,1)
pvec <- seq(0.1, 0.9, 0.2)
tolerance <- 0.01

test_that("pml() and qml() are compatible for Type 1", {
  for (p in pvec){
    for (scale in scalevec){
      for (alpha in alphavec){
        q <- qml(p = p, alpha = alpha, scale = scale, second.type = FALSE)
        p2 <- pml(q = q, alpha = alpha, scale = scale, second.type = FALSE)
        expect_equal(object = p2, expected = p, tolerance = tolerance)
      }
    }
  }
})

test_that("pml() and qml() are compatible for Type 2", {
  for (p in pvec){
    for (scale in scalevec){
      for (alpha in alphavec){
        q <- qml(p = p, alpha = alpha, scale = scale, second.type = TRUE)
        p2 <- pml(q = q, alpha = alpha, scale = scale, second.type = TRUE)
        expect_equal(object = p2, expected = p, tolerance = tolerance)
      }
    }
  }
})