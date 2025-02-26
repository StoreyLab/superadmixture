## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----status, echo=FALSE-------------------------------------------------------
fast_run <- FALSE

## ----load_amr-----------------------------------------------------------------
data("X_amr",   package = "superadmixture")
data("fam_amr", package = "superadmixture")

## ----estimate_coanc_indiv_amr, eval=!fast_run, message=FALSE, warning=FALSE----
# load popkin package
library(popkin)

# estimate kinship coefficients
# `popkin` function requires a m (number of loci) x n (number of individuals) as input, so we transpose the `X_amr` first 
kinship <- popkin(t(X_amr)) 

# map kinship coefficients to coancestry coefficients
coanc_indiv <- inbr_diag(kinship)

# kinship coefficients and coancestry coefficients are probability, 
# so we truncate at zero to avoid negative probability
coanc_indiv <- ifelse(coanc_indiv < 0, 0, coanc_indiv)

## ----helper_functions, eval=TRUE, warning=FALSE, message=FALSE----------------
plot_colors_subpops <- function(colors, y = FALSE) {
  # number of individuals
  n <- length(colors)
  if (y) {
    x <- rbind(n: 1)
  } else {
    x <- cbind(1: n)
  }
  image(x, col = colors, axes = FALSE, useRaster = TRUE)
}

# We also need to construct a legend for these colors.
legend_color_categories <- function(colors, categories, label, cex_label = 1) {
  x <- 1: length(colors)
  image(y = x, z = rbind(x), col = colors, xaxt = "n", yaxt = "n")
  axis(4, at = x, labels = categories, tick = FALSE)
  mtext(side = 4, label, line = 2, cex = cex_label)
}

## ----plot_coanc_indiv_os_amr, eval=TRUE, warning=FALSE, message=FALSE, fig.width=6, fig.height=5.3, fig.align='center'----
# We plot the individual-level coancestry 
subpop_order <- c('PUR', 'CLM', 'PEL', 'MXL')
n_subpops      <- length(subpop_order)
colors_subpops <- RColorBrewer::brewer.pal(n_subpops, "Set3")
fam_amr$col    <- colors_subpops[match(fam_amr$fam, subpop_order)]

par(mar = c(0, 0, 1, 0) + 0.2)

layout(rbind(c(3, 1, 2), c(3, 1, 5), c(0, 4, 0)), widths  = c(0.05,   1,  0.2), heights = c(0.5,  0.5, 0.07))

popkin::plot_popkin(kinship = coanc_indiv, layout_add = FALSE, ylab = '', leg_title = "Coancestry", titles = c("OS individual coancestry"))

plot_colors_subpops(fam_amr$col, y = TRUE)

mtext('Individuals', side = 2, line = 0.5, xpd = NA, cex = 0.8)

plot_colors_subpops(fam_amr$col)

par(mar = c(0, 0, 0, 3) + 0.2)
legend_color_categories(colors = colors_subpops, categories = subpop_order, label = 'Subpopulation', cex_label = 0.8)

## ----estimate_admix_props_amr, eval=!fast_run, message=FALSE, warning=FALSE----
library(superadmixture)

# estimate individual-specific allele frequencies
obj <- est_p_indiv(X_amr, k_antepops = 3, loci_on_cols = TRUE)
p_indiv  <- obj$p_indiv
rowspace <- obj$rowspace

# estimate P and Q by decomposing individual-specific allele frequencies
obj <- factor_p_indiv(p_indiv, k_antepops = 3, rowspace = rowspace, verbose = FALSE, max_iters = 200, tol = 1e-2)
Q_hat <- obj$Q_hat

## ----estimate_coanc_pops_amr, eval=!fast_run, message=FALSE, warning=FALSE----
# estimate population coancestry under the super admixture model
coanc_pops_sup <- est_coanc(coanc_indiv, Q_hat, model = "super")

# estimate population coancestry under the standard admixture model
coanc_pops_std <- est_coanc(coanc_indiv, Q_hat, model = "standard")

## ----fit_tree_amr, eval=TRUE, message=FALSE, warning=FALSE--------------------
# reorder antecedent populations in the ascending order of coancestry
index <- order(diag(coanc_pops_sup))
Q_hat <- Q_hat[index, ]
coanc_pops_sup <- coanc_pops_sup[index, index]

# label antecedent populations
colnames(coanc_pops_sup) <- rownames(coanc_pops_sup) <- paste0("S", 1:3)

# fit tree using the `fit_tree` function from the `bnpsd` package
tree <- bnpsd::fit_tree(coanc_pops_sup)

## ----draw_tree_amr, eval=FALSE, warning=FALSE, message=FALSE, fig.width=8, fig.height=3, fig.align="center"----
#  colors <- get_seq_colors("Reds", 3)
#  names(colors) <- paste0("S", 1:3)
#  fig_tree <- plot_tree(tree, colors = colors, font_size = 17)
#  fig_tree

## ----load_tree_amr, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE-------
library(knitr)
include_graphics("figures/tree-amr.png")

## ----draw_admix_props_amr, eval=TRUE, warning=FALSE, message=FALSE, fig.width=8, fig.height=3, fig.align="center"----
colors <- get_seq_colors("Reds", 3)
names(colors) <- paste0("S", 1:3)
barplot_admix(Q_hat, colors = colors, indiv_on_cols = TRUE)

## ----draw_coanc_antepops_amr, eval=TRUE, warning=FALSE, message=FALSE, fig.width=8, fig.height=3, fig.align="center"----
par(xpd = TRUE)
heatmap_coanc_antepops(coanc_pops_sup, tl.offset = 1)

## ----estimate_coanc_supadmix_stdadmix_amr, eval=!fast_run, message=FALSE, warning=FALSE----
coanc_sup <- t(Q_hat) %*% coanc_pops_sup %*% Q_hat
coanc_std <- t(Q_hat) %*% coanc_pops_std %*% Q_hat

## ----plot_coanc_indiv_sup_amr, eval=TRUE, warning=FALSE, message=FALSE, fig.width=6, fig.height=5.3, fig.align='center'----
par(mar = c(0, 0, 1, 0) + 0.2)

layout(rbind(c(3, 1, 2), c(3, 1, 5), c(0, 4, 0)), widths  = c(0.05,   1,  0.2), heights = c(0.5,  0.5, 0.07))

popkin::plot_popkin(kinship = coanc_sup, layout_add = FALSE, leg_cex = 0.8, ylab = '', leg_title = "Coancestry", titles = c("Super admixture individual coancestry"))

plot_colors_subpops(fam_amr$col, y = TRUE)

mtext('Individuals', side = 2, line = 0.5, xpd = NA, cex = 0.8)

plot_colors_subpops(fam_amr$col)

par(mar = c(0, 0, 0, 3) + 0.2)
legend_color_categories(colors = colors_subpops, categories = subpop_order, label = 'Subpopulation', cex_label = 0.8)

## ----plot_coanc_indiv_std_amr, eval=TRUE, warning=FALSE, message=FALSE, fig.width=6, fig.height=5.3, fig.align='center'----
par(mar = c(0, 0, 1, 0) + 0.2)

layout(rbind(c(3, 1, 2), c(3, 1, 5), c(0, 4, 0)), widths  = c(0.05,   1,  0.2), heights = c(0.5,  0.5, 0.07))

popkin::plot_popkin(kinship = coanc_std, layout_add = FALSE, leg_cex = 0.8, ylab = '', leg_title = "Coancestry", titles = c("Standard admixture individual coancestry"))

plot_colors_subpops(fam_amr$col, y = TRUE)

mtext('Individuals', side = 2, line = 0.5, xpd = NA, cex = 0.8)

plot_colors_subpops(fam_amr$col)

par(mar = c(0, 0, 0, 3) + 0.2)
legend_color_categories(colors = colors_subpops, categories = subpop_order, label = 'Subpopulation', cex_label = 0.8)

## ----sim_geno_amr, eval=!fast_run, message=FALSE, warning=FALSE---------------
# approximate ancestral allele frequencies by average allele frequencies
p_anc <- 0.5 * colMeans(X_amr, na.rm = TRUE)
# simulate genotypes according to the double-admixture method
X_sim_amr <- dbl_admixture(p_anc, coanc_pops_sup, Q_hat, geno_only = TRUE)

## ----estimate_sim_coanc_indiv_amr, eval=!fast_run, message=FALSE, warning=FALSE----
# estimate the kinship of simulated genotype
# since OS method assumes minimum of pairwise kinship, which doesn't hold 
# we use the following strategy to adjust OS kinship estimate
kinship_sim  <- popkin(X_sim_amr)
kinship_sim  <- (kinship_sim - 1) * (1 - min(coanc_sup[col(coanc_sup) != row(coanc_sup)])) + 1

# map kinship coefficients to coancestry coefficients
coanc_sim_indiv   <- inbr_diag(kinship_sim)

# kinship coefficients and coancestry coefficients are probability, 
# so we truncate at zero to avoid negative probability
kinship_sim <- ifelse(kinship_sim < 0, 0, kinship_sim)
coanc_sim_indiv  <- ifelse(coanc_sim_indiv < 0, 0, coanc_sim_indiv)

## ----plot_coanc_indiv_os_sim_amr, eval=TRUE, warning=FALSE, message=FALSE, fig.width=6, fig.height=5.3, fig.align='center'----
par(mar = c(0, 0, 1, 0) + 0.2)

layout(rbind(c(3, 1, 2), c(3, 1, 5), c(0, 4, 0)), widths  = c(0.05,   1,  0.2), heights = c(0.5,  0.5, 0.07))

popkin::plot_popkin(kinship = coanc_sim_indiv, layout_add = FALSE, leg_cex = 0.8, ylab = '', leg_title = "Coancestry", titles = c("OS individual coancestry of simulated data"))

plot_colors_subpops(fam_amr$col, y = TRUE)

mtext('Individuals', side = 2, line = 0.5, xpd = NA, cex = 0.8)

plot_colors_subpops(fam_amr$col)

par(mar = c(0, 0, 0, 3) + 0.2)
legend_color_categories(colors = colors_subpops, categories = subpop_order, label = 'Subpopulation', cex_label = 0.8)

## ----load_hgdp----------------------------------------------------------------
data("X_hgdp",   package = "superadmixture")
data("fam_hgdp", package = "superadmixture")

## ----estimate_coanc_indiv_hgdp, eval=!fast_run, message=FALSE, warning=FALSE----
# we provide the population labels
kinship <- popkin(t(X_hgdp), subpops = fam_hgdp$population) 
coanc_indiv <- inbr_diag(kinship)
kinship <- ifelse(kinship< 0, 0, kinship)
coanc_indiv <- ifelse(coanc_indiv < 0, 0, coanc_indiv)

## -----------------------------------------------------------------------------
plot_colors_subpops <- function(pops, srt = 0, cex = 0.6, y = FALSE) {
  n <- length(pops)
  k <- unique(pops)
  pops <- factor(pops, levels = unique(pops))
  xintercept <- cumsum(table(pops))
  breaks <- xintercept - 0.5 * as.numeric(table(pops))
  if (y) {
    plot(NULL, xlim = c(0, 1), ylim = c(1, n), axes = FALSE, ann = FALSE, xaxs = "i", yaxs = "i")
    text(1, n-rev(breaks), rev(unique(pops)), cex = cex, srt = srt, xpd = TRUE, adj = c(1, 0.5))
  } else {
    plot(NULL, xlim = c(1, n), ylim = c(0, 1), axes = FALSE, ann = FALSE, xaxs = "i", yaxs = "i")
    text(breaks, 1, unique(pops), cex = cex, srt = srt, xpd = TRUE, adj = c(1, 0.5))
  }
}

## ----plot_coanc_indiv_os_hgdp, eval=TRUE, warning=FALSE, message=FALSE, fig.width=7, fig.height=5.3, fig.align='center'----
par(mar = c(0, 0, 0, 0) + 0.2)
layout(rbind(c(3, 1, 2), c(3, 1, 5), c(0, 4, 0)), widths = c(0.1, 1, 0.2), heights = c(0.5, 0.5, 0.1))

popkin::plot_popkin(kinship = coanc_indiv, layout_add = FALSE, leg_cex = 0.8, labs_text = FALSE, labs_lwd = 0.1, labs = fam_hgdp$subpop, ylab = '', leg_title = "Coancestry")

par(mar = c(0.2, 0, 0.2, 0))
plot_colors_subpops(fam_hgdp$subpop, y = TRUE, cex = 0.8)
mtext('Individuals', side = 2, line = 0.5, xpd = NA, cex = 0.8)
par(mar = c(0, 0.2, 0, 0.2))
plot_colors_subpops(fam_hgdp$subpop, srt = 90, cex = 0.8)

## ----estimate_admix_props_hgdp, eval=!fast_run, message=FALSE, warning=FALSE----
library(superadmixture)

# estimate individual-specific allele frequencies
obj <- est_p_indiv(X_hgdp, k_antepops = 7, loci_on_cols = TRUE)
p_indiv <- obj$p_indiv
rowspace <- obj$rowspace

# estimate P and Q by decomposing individual-specific allele frequencies
obj <- factor_p_indiv(p_indiv, k_antepops = 7, rowspace = rowspace, verbose = FALSE, max_iters = 200, tol = 1e-2)
Q_hat <- obj$Q_hat

## ----estimate_coanc_antepops_hgdp, eval=!fast_run, message=FALSE, warning=FALSE----
# estimate population coancestry under the super admixture model
coanc_pops_sup <- est_coanc(coanc_indiv, Q_hat, model = "super")

# estimate population coancestry under the standard admixture model
coanc_pops_std <- est_coanc(coanc_indiv, Q_hat, model = "standard")

## ----fit_tree_hgdp, eval=TRUE, message=FALSE, warning=FALSE-------------------
# reorder antecedent populations in the ascending order of coancestry
index <- order(diag(coanc_pops_sup))
Q_hat <- Q_hat[index, ]
coanc_pops_sup <- coanc_pops_sup[index, index]

# label antecedent populations
colnames(coanc_pops_sup) <- rownames(coanc_pops_sup) <- paste0("S", 1:7)

# fit tree using the `fit_tree` function from the `bnpsd` package
tree <- bnpsd::fit_tree(coanc_pops_sup)

## ----draw_uncolorred_tree_hgdp, eval=FALSE, warning=FALSE, message=FALSE, fig.width=8, fig.height=3, fig.align="center"----
#  fig_tree <- plot_tree(tree, font_size = 17)
#  fig_tree

## ----load_tree_uncolorred_hgdp, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE----
include_graphics("figures/tree-uncolorred-hgdp.png")

## ----draw_tree_hgdp, eval=FALSE, warning=FALSE, message=FALSE, fig.width=8, fig.height=3, fig.align="center"----
#  colors <- c(get_seq_colors("Blues", 2), get_seq_colors("Greens", 2), get_seq_colors("Reds", 2), get_seq_colors("Purples", 1))
#  names(colors) <- paste0("S", 1:7)
#  fig_tree <- plot_tree(tree, colors = colors, font_size = 17)
#  fig_tree

## ----load_tree_hgdp, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE------
include_graphics("figures/tree-hgdp.png")

## ----draw_admix_props_hgdp, eval=TRUE, warning=FALSE, message=FALSE, fig.width=8, fig.height=3, fig.align="center"----
colors <- c(get_seq_colors("Blues", 2), get_seq_colors("Greens", 2), get_seq_colors("Reds", 2), get_seq_colors("Purples", 1))
names(colors) <- paste0("S", 1:7)
barplot_admix(Q_hat, colors = colors, subpops = fam_hgdp$subpop, indiv_on_cols = TRUE)

## ----draw_coanc_antepops_hgdp, eval=TRUE, warning=FALSE, message=FALSE, fig.width=8, fig.height=3, fig.align="center"----
par(xpd = TRUE)
heatmap_coanc_antepops(coanc_pops_sup, tl.offset = 1)

## ----estimate_coanc_supadmix_stdadmix_hgdp, eval=!fast_run, message=FALSE, warning=FALSE----
coanc_sup <- t(Q_hat) %*% coanc_pops_sup %*% Q_hat
coanc_std <- t(Q_hat) %*% coanc_pops_std %*% Q_hat

## ----plot_coanc_indiv_sup_hgdp, eval=TRUE, warning=FALSE, message=FALSE, fig.width=7, fig.height=5.3, fig.align='center'----
par(mar = c(0, 0, 0, 0) + 0.2)
layout(rbind(c(3, 1, 2), c(3, 1, 5), c(0, 4, 0)), widths = c(0.1, 1, 0.2), heights = c(0.5, 0.5, 0.1))

popkin::plot_popkin(kinship = coanc_sup, layout_add = FALSE, leg_cex = 0.8, labs_text = FALSE, labs_lwd = 0.1, labs = fam_hgdp$subpop, ylab = '', leg_title = "Coancestry")

par(mar = c(0.2, 0, 0.2, 0))
plot_colors_subpops(fam_hgdp$subpop, y = TRUE, cex = 0.8)
mtext('Individuals', side = 2, line = 0.5, xpd = NA, cex = 0.8)
par(mar = c(0, 0.2, 0, 0.2))
plot_colors_subpops(fam_hgdp$subpop, srt = 90, cex = 0.8)

## ----plot_coanc_indiv_std_hgdp, eval=TRUE, warning=FALSE, message=FALSE, fig.width=7, fig.height=5.3, fig.align='center'----
par(mar = c(0, 0, 0, 0) + 0.2)
layout(rbind(c(3, 1, 2), c(3, 1, 5), c(0, 4, 0)), widths = c(0.1, 1, 0.2), heights = c(0.5, 0.5, 0.1))

popkin::plot_popkin(kinship = coanc_std, layout_add = FALSE, leg_cex = 0.8, labs_text = FALSE, labs_lwd = 0.1, labs = fam_hgdp$subpop, ylab = '', leg_title = "Coancestry")

par(mar = c(0.2, 0, 0.2, 0))
plot_colors_subpops(fam_hgdp$subpop, y = TRUE, cex = 0.8)
mtext('Individuals', side = 2, line = 0.5, xpd = NA, cex = 0.8)
par(mar = c(0, 0.2, 0, 0.2))
plot_colors_subpops(fam_hgdp$subpop, srt = 90, cex = 0.8)

## ----sim_geno_hgdp, eval=!fast_run, message=FALSE, warning=FALSE--------------
# approximate ancestral allele frequencies by average allele frequencies
p_anc <- 0.5 * colMeans(X_hgdp, na.rm = TRUE)
# simulate genotypes according to the double-admixture method
X_sim_hgdp <- dbl_admixture(p_anc, coanc_pops_sup, Q_hat, geno_only = TRUE)

## ----estimate_sim_coanc_indiv_hgdp, eval=!fast_run, message=FALSE, warning=FALSE----
# estimate the kinship of simulated genotype
# since OS method assumes minimum of pairwise kinship, which doesn't hold 
# we use the following strategy to adjust OS kinship estimate
kinship_sim  <- popkin(X_sim_hgdp)
kinship_sim  <- (kinship_sim - 1) * (1 - min(coanc_sup[col(coanc_sup) != row(coanc_sup)])) + 1

# map kinship coefficients to coancestry coefficients
coanc_sim_indiv   <- inbr_diag(kinship_sim)

# kinship coefficients and coancestry coefficients are probability, 
# so we truncate at zero to avoid negative probability
kinship_sim <- ifelse(kinship_sim < 0, 0, kinship_sim)
coanc_sim_indiv   <- ifelse(coanc_sim_indiv < 0, 0, coanc_sim_indiv)

## ----plot_coanc_indiv_sim_os_hgdp, eval=TRUE, warning=FALSE, message=FALSE, fig.width=7, fig.height=5.3, fig.align='center'----
par(mar = c(0, 0, 0, 0) + 0.2)
layout(rbind(c(3, 1, 2), c(3, 1, 5), c(0, 4, 0)), widths = c(0.1, 1, 0.2), heights = c(0.5, 0.5, 0.1))

popkin::plot_popkin(coanc_sim_indiv, layout_add = FALSE, leg_cex = 0.8, labs_text = FALSE, labs_lwd = 0.1, labs = fam_hgdp$subpop, ylab = '', leg_title = "Coancestry")

par(mar = c(0.2, 0, 0.2, 0))
plot_colors_subpops(fam_hgdp$subpop, y = TRUE, cex = 0.8)
mtext('Individuals', side = 2, line = 0.5, xpd = NA, cex = 0.8)
par(mar = c(0, 0.2, 0, 0.2))
plot_colors_subpops(fam_hgdp$subpop, srt = 90, cex = 0.8)

