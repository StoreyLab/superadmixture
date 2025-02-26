test_that("rdirichlet() errors if the input parameter(s) are not valid", {
  expect_error(rdirichlet(2, c(0.5, 0)))
  expect_error(rdirichlet(0, c(0.5, 0.5)))
  expect_error(rdirichlet(3, matrix(1:4, 2, 2)))
})

test_that("rdirichlet() can generate samples from Dirichlet distribution", {
  EPSILON <- 1e-2
  n <- 1000000
  alpha <- runif(3)
  samples <- rdirichlet(n, alpha)
  expect_lt(norm(colMeans(samples) - alpha / sum(alpha), "2"), EPSILON)
  expect_lt(norm(apply(samples, 2, var) - alpha * (sum(alpha) - alpha) / (sum(alpha)^2 * (sum(alpha) + 1)), "2"), EPSILON)

  alpha <- sample(1:10, 5)
  samples <- rdirichlet(n, alpha)
  expect_lt(norm(colMeans(samples) - alpha / sum(alpha), "2"), EPSILON)
  expect_lt(norm(apply(samples, 2, var) - alpha * (sum(alpha) - alpha) / (sum(alpha)^2 * (sum(alpha) + 1)), "2"), EPSILON)
})
