library(BEDMatrix)
library(genio)
library(devtools)
library(popkin)
library(seriation)
library(superadmixture)

set.seed(525)

X_amr    <- BEDMatrix::BEDMatrix("amr-rawdata/amr_final", simple_names = TRUE)
fam_amr  <- genio::read_fam( "amr-rawdata/amr_final")

X_amr  <- X_amr[, sample(1:ncol(X_amr), size = 2000, replace = FALSE)]

# Re-order individuals by population label
subpop_order <- c('PUR', 'CLM', 'PEL', 'MXL')
index        <- order(factor(fam_amr$fam, level = subpop_order))
fam_amr      <- fam_amr[index, ]
X_amr       <- X_amr[index, ]

# Estimate the Ochoa-Storey kinship
obj         <- popkin::popkin_A(t(X_amr))
A_min       <- popkin::popkin_A_min_subpops(obj$A)
kinship_os  <- 1 - obj$A / A_min

# Order individuals using seriation
order_seriate <- function(kinship) {
  distance <- -kinship
  distance <- distance - min(distance)
  distance <- as.dist(distance)

  # perform desired optimization
  seriation_object <- seriation::seriate(distance, method = 'ARSA')
  index <- seriation::get_order(seriation_object)
  y <- diag(kinship)[index]
  x <- 1:length(y)
  m <- coef(lm(y ~ x))[2]
  if (m < 0) {
    index <- rev(index)
  }
  return(index)
}

index      <- order_seriate(kinship_os)
fam_amr    <- fam_amr[index, ]
X_amr      <- X_amr[index, ]

# estimate kinship coefficients
kinship <- popkin(t(X_amr))
coanc_indiv <- inbr_diag(kinship)
coanc_indiv <- ifelse(coanc_indiv < 0, 0, coanc_indiv)

# estimate admixture proportions
obj <- est_p_indiv(X_amr, k_antepops = 3, loci_on_cols = TRUE)
p_indiv  <- obj$p_indiv
rowspace <- obj$rowspace
obj <- factor_p_indiv(p_indiv, k_antepops = 3, rowspace = rowspace, verbose = FALSE, max_iters = 200, tol = 1e-2)
admix_props_amr <- obj$Q_hat

# estimate population coancestry under the super-admixture model
coanc_pops_amr <- est_coanc(coanc_indiv, admix_props_amr, model = "super")

# approximate ancestral allele frequencies by average allele frequencies
p_anc_amr <- 0.5 * colMeans(X_amr, na.rm = TRUE)

use_data(X_amr, fam_amr, admix_props_amr, coanc_pops_amr, p_anc_amr, overwrite = TRUE)
