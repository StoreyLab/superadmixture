test_that("check_p_anc() errors if the input is not numeric", {
  expect_error(check_p_anc(c("A", "B", "C", "D")))
})

test_that("check_p_anc() errors if the input is not a vector", {
  expect_error(check_p_anc(matrix(1:4, 2, 2)))
})

test_that("check_p_anc() errors if some entries are greater than 1", {
  expect_error(check_p_anc(c(1.5, 0.4, 0.3)))
})

test_that("check_p_anc() errors if some entries are less than 0", {
  expect_error(check_p_anc(c(-0.5, 0.4, 0.3)))
})

test_that("check_coancestry() errors if the input is not symmetric", {
  expect_error(check_coancestry(matrix(runif(10), 5, 2)))
  expect_error(check_coancestry(matrix(runif(10), 2, 5)))
})

test_that("check_coancestry() errors if the input is not numeric", {
  expect_error(check_coancestry(matrix(c("A", "B", "C", "D"), 2, 2)))
})

test_that("check_coancestry() works if the input is a valid coancestry", {
  k_antepops <- sample(3:9, 1)
  expect_no_condition(check_coancestry(construct_coanc(k_antepops, "1")))
  expect_no_condition(check_coancestry(construct_coanc(k_antepops, "2")))
  expect_no_condition(check_coancestry(construct_coanc(k_antepops, "3")))
})

test_that("check_coancestry() errors if some entries are less than 0", {
  A <- matrix(c(1, -1, -1, 1), 2, 2)
  expect_error(check_coancestry(A))
})

test_that("check_coancestry() errors if some entries are greater than 1", {
  A <- matrix(c(1.5, 0, 0, 1.5), 2, 2)
  expect_error(check_coancestry(A))
})

test_that("check_admix_proportions() errors if the input is not numeric", {
  expect_error(check_admix_proportions(matrix(c("A", "B", "C", "D"), 2, 2)))
})

test_that("check_admix_proportions() errors if some of the entries are negative", {
  A <- matrix(c(1, 1, -1, -1, 1, 1), 3, 2)
  expect_error(check_admix_proportions(A))
})

test_that("check_admix_proportions() errors if `rowSums(admix_proportions) != 1` && `colSums(admix_proportions) != 1`", {
  A <- matrix(c(1, 2, 3, 4), 2, 2)
  expect_error(check_admix_proportions(A))
})

test_that("check_admix_proportions() works if the input is valid", {
  k_antepops <- 3
  expect_no_condition(check_admix_proportions(construct_admix_props(alpha_type = "1")))
  expect_no_condition(check_admix_proportions(construct_admix_props(alpha_type = "2")))
  expect_no_condition(check_admix_proportions(construct_admix_props(alpha_type = "3")))
  expect_no_condition(check_admix_proportions(construct_admix_props(alpha_type = "spatial")))
})
