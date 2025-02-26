# check if `p_anc` is valid
check_p_anc <- function(p_anc, EPSILON = 1e-3) {
  stopifnot("`p_anc` is a numeric vector" = is.numeric(p_anc) && is.vector(p_anc))
  stopifnot("all entries of `p_anc` are less than 1.0" = all(p_anc - 1 < EPSILON))
  stopifnot("all entries of `p_anc` are greater than 0.0" = all(p_anc > -EPSILON))
}

# check if `p_antepops` is valid
check_p_antepops <- function(p_antepops, EPSILON = 1e-3) {
  stopifnot("`p_antepops` is a numeric matrix" = is.numeric(p_antepops) && is.vector(p_antepops))
  stopifnot("all entries of `p_antepops` are less than 1.0" = all(p_antepops - 1 < EPSILON))
  stopifnot("all entries of `p_antepops` are greater than 0.0" = all(p_antepops > -EPSILON))
}

# check if `p_indiv` is valid
check_p_indiv <- function(p_indiv, EPSILON = 1e-3) {
  stopifnot("`p_indiv` is a numeric matrix" = is.numeric(p_indiv) && is.vector(p_indiv))
  stopifnot("all entries of `p_indiv` are less than 1.0" = all(p_indiv - 1 < EPSILON))
  stopifnot("all entries of `p_indiv` are greater than 0.0" = all(p_indiv > -EPSILON))
}

# check if `coanc_individuals` is valid
check_coancestry <- function(coancestry, EPSILON = 1e-3) {
  stopifnot("`coancestry` is a numeric matrix" = is.numeric(coancestry) && is.matrix(coancestry))
  stopifnot("`coancestry` is a symmetric matrix" = isSymmetric(coancestry))
  stopifnot("All entries of `coancestry` are less than 1.0" = all(coancestry - 1 < EPSILON))
  stopifnot("All entries of `coancestry` are greater than 0" = all(coancestry > -EPSILON))
}

# check if `coanc_individuals` is valid
check_admix_proportions <- function(admix_proportions, EPSILON = 1e-3) {
  stopifnot("`admix_proportions` is a numeric matrix" = is.numeric(admix_proportions) && is.matrix(admix_proportions))
  stopifnot("some entries of `admix_proportions` are greater than 1.0" = all(admix_proportions > -EPSILON))
  stopifnot(
    "`rowSums(admix_proportions) = 1` || `colSums(admix_proportions) = 1`" =
      all(abs(rowSums(admix_proportions) - 1) < EPSILON) || all(abs(colSums(admix_proportions) - 1) < EPSILON)
  )
}

# check if `X` is valid
check_geno <- function(X) {
  stopifnot("`X` takes values in 0, 1, 2, NA" = all(X %in% c(0, 1, 2, NA)))
}
