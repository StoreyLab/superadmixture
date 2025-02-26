library(BEDMatrix)
library(genio)
library(devtools)
library(dplyr)
library(readr)
library(popkin)
library(superadmixture)

set.seed(525)

X_hgdp    <- BEDMatrix("hgdp-rawdata/hgdp_wgs_autosomes_final", simple_names = TRUE)
fam_hgdp  <- read_fam( "hgdp-rawdata/hgdp_wgs_autosomes_final")
info_hgdp <- read_tsv( "hgdp-rawdata/hgdp_wgs.20190516.metadata.txt", col_types = 'ccccccddccddddd') %>%
  select(c('sample', 'population', 'latitude', 'longitude', 'region')) 

fam_hgdp <- left_join(fam_hgdp, info_hgdp, by = c("id" = "sample")) %>%
            rename("subpop" = "region") %>%
            mutate(subpop = ifelse(subpop == "AFRICA", "Africa",
                            ifelse(subpop == "MIDDLE_EAST", "MiddleEast",
                            ifelse(subpop == "EUROPE", "Europe",
                            ifelse(subpop == "CENTRAL_SOUTH_ASIA", "CSAsia",
                            ifelse(subpop == "EAST_ASIA", "EAsia",
                            ifelse(subpop == "AMERICA", "Americas", "Oceania")))))))


# reorder by subpop
subpop_order <- c('Africa', 'MiddleEast', 'Europe', 'CSAsia', 'EAsia', 'Americas', 'Oceania')
index <- order( match(fam_hgdp$subpop, subpop_order) )
fam_hgdp   <- fam_hgdp[index, ]
X_hgdp     <- X_hgdp[index, ]

# subset loci
X_hgdp  <- X_hgdp[, sample(1:ncol(X_hgdp), size = 10000, replace = FALSE)]

# estimate coancestry among individuals
kinship <- popkin(t(X_hgdp), subpops = fam_hgdp$population) 
coanc_indiv <- inbr_diag(kinship)
coanc_indiv <- ifelse(coanc_indiv < 0, 0, coanc_indiv)

# estimate individual-specific allele frequencies
obj <- est_p_indiv(X_hgdp, k_antepops = 7, loci_on_cols = TRUE)
p_indiv <- obj$p_indiv
rowspace <- obj$rowspace

# estimate P and Q by decomposing individual-specific allele frequencies
obj <- factor_p_indiv(p_indiv, k_antepops = 7, rowspace = rowspace, verbose = FALSE, max_iters = 200, tol = 1e-2)
admix_props_hgdp <- obj$Q_hat

# estimate population coancestry under the super-admixture model
coanc_pops_hgdp <- est_coanc(coanc_indiv, admix_props_hgdp, model = "super")

# approximate ancestral allele frequencies by average allele frequencies
p_anc_hgdp <- 0.5 * colMeans(X_hgdp, na.rm = TRUE)

use_data(X_hgdp, fam_hgdp, admix_props_hgdp, coanc_pops_hgdp, p_anc_hgdp, overwrite = TRUE)
