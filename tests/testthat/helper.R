construct_coanc <- function(k_antepops = 3, coanc_type = "1") {
  switch(coanc_type,
         "1" = {
           A <- matrix(runif(k_antepops^2, min = 0, max = 0.3), ncol = k_antepops)
           return(t(A) %*% A)
         },
         "2" = {
           Gamma <- diag(runif(k_antepops, min = 0, max = 1))
           W <- t(rdirichlet(n = k_antepops, a = rep(1, k_antepops)))
           return(t(W) %*% Gamma %*% W)
         },
         "3" = {
           Gamma <- diag(runif(10 * k_antepops, min = 0, max = 1))
           W <- t(rdirichlet(n = k_antepops, a = rep(1, 10 * k_antepops)))
           return(t(W) %*% Gamma %*% W)
         }
  )
}

construct_admix_props <- function(k_antepops = 3, n_inds = 1000, alpha_type = "1") {
  alpha <- switch(alpha_type,
    "1" = {
      prototype <- c(1, 1, 1)
      c(rep(prototype, each = k_antepops %/% length(prototype)),
        rep(prototype, length.out = k_antepops %% length(prototype)))
    },
    "2" = {
      prototype <- c(.1, .1, .1)
      c(rep(prototype, each = k_antepops %/% length(prototype)),
        rep(prototype, length.out = k_antepops %% length(prototype)))
    },
    "3" = {
      prototype <- c(2, 1, .5)
      c(rep(prototype, each = k_antepops %/% length(prototype)),
        rep(prototype, length.out = k_antepops %% length(prototype)))
    },
    "spatial" = {"spatial"}
  )
  if (is.numeric(alpha)) {
    return(t(rdirichlet(n = n_inds, alpha = alpha)))
  } else {
    return(t(bnpsd::admix_prop_1d_linear(n_inds, k_antepops, sigma = 0.5)))
  }
}


D_binomial <- function(X) {
  s <- 2
  m <- dim(X)[1]
  n <- dim(X)[2]
  delta_hat <- rep(0, n)
  v <- function(x, s) (s * x - x^2)/(s - 1)
  for (i in 1:n) {
    delta_hat[i] <- 1/m * sum(v(X[, i], s))
  }
  D <- diag(delta_hat)
}

lse <- function (X, d, svd_method = "base"){
  m <- dim(X)[1]
  D <- D_binomial(X)
  if (svd_method == "base") {
    rowspace <- eigen(1/m * t(X) %*% X - D)
    vectors <- rowspace$vectors[, 1:d]
    values <- rowspace$values[1:d]
  }
  else if (svd_method == "truncated_svd") {
    rowspace <- svd::propack.svd(1/m * t(X) %*% X - D, neig = d)
    vectors <- rowspace$v[, 1:d]
    values <- rowspace$d[1:d]
  }
  vals <- list(vectors = vectors, values = values)
  return(vals)
}

estimate_F <- function(X, d, svd_method = "base") {
    m <- dim(X)[1]
    D <- D_binomial(X)
    rowspace <- lse(X = X, d = d, svd_method = svd_method)
    F_hat <- 1/2 * X %*% rowspace$vectors[, 1:d] %*% t(rowspace$vectors[,1:d])
    F_hat[F_hat < 0] <- 0
    F_hat[F_hat > 1] <- 1
    vals <- list(F_hat = F_hat, rowspace = rowspace)
}
