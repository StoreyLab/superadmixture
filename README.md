# The super admixture framework

![R-CMD-check.yaml](https://github.com/StoreyLab/superadmixture/actions/workflows/R-CMD-check.yaml/badge.svg)

The `superadmixture` R package implements functions for inferring the population coancestry and simulating genotypes according to the super admixture framework.  It also has functions for conducting the structured Hardy-Weinberg equilibrium test and estimating admixture proportions using the ALStructure algorithm. 

## Installation

```R
install.packages("devtools")
devtools::install_github("StoreyLab/superadmixture")
```

## Examples 

### Input

We created a subset of HGDP datasets by first applying allele frequency filters and then LD-pruning to the HGDP dataset. We then randomly selected 10,000 SNPs out of LD-pruned SNP sets. This subset of genotypes is available in `data/X_hgdp.rda`. This data has 929 individuals and 10,000 loci. The associated fam file can be found in `data/fam_hgdp.rda`. These data can be reproduced by scripts `data-raw/{hgdp.bash,hgdp.R}`. 

```R
data("X_hgdp",   package = "superadmixture")
data("fam_hgdp", package = "superadmixture")
```

### Estimating oancestry among individuals

```R
library(popkin)

# estimate kinship coefficients
kinship <- popkin(t(X_hgdp)) 

# map kinship coefficients to coancestry coefficients
coanc_indiv  <- inbr_diag(kinship)

# kinship coefficients and coancestry coefficients are probability, 
# so we truncate at zero to avoid negative probability
coanc_indiv   <- ifelse(coanc_indiv < 0, 0, coanc_indiv)
```

### Estimating admixture proportions

```R
library(superadmixture)

# estimate individual-specific allele frequencies
obj <- est_p_indiv(X_hgdp, k_antepops = 7, loci_on_cols = TRUE)
p_indiv  <- obj$p_indiv
rowspace <- obj$rowspace

# estimate P and Q by decomposing individual-specific allele frequencies
obj <- factor_p_indiv(p_indiv, k_antepops = 7, rowspace = rowspace, verbose = FALSE, max_iters = 200, tol = 1e-2)
Q_hat <- obj$Q_hat
```

### Estimating coancestry among antecedent populations 

After obtaining `coanc_indiv` and `Q_hat`, we can use the function `est_coanc` to estimate population coancestry under the super admixture model and under the standard admixture model. 

```R
# estimate population coancestry under the super-admixture model
coanc_pops_sup <- est_coanc(coanc_indiv, Q_hat, model = "super")

# estimate population coancestry under the standard admixture model
coanc_pops_std <- est_coanc(coanc_indiv, Q_hat, model = "standard")
``` 

### Visualizing population-level coancestry and admixture proportions

```R
# order `Q_hat` and `coanc_pops_sup` in the ascending order of coancestry
index <- order(diag(coanc_pops_sup))
coanc_pops_sup <- coanc_pops_sup[index, index]
Q_hat <- Q_hat[index, ]
```

```R
# use the `fit_tree` function from bnpsd package to create a phylogenetic tree
library(bnpsd)
colnames(coanc_pops_sup) <- rownames(coanc_pops_sup) <- paste0("S", 1:7)
tree <- fit_tree(coanc_pops_sup)
```

We can visualize `tree` using the `plot_tree` function in the superadmixture package. Based on the topology of the tree (see the package vignette for details), we decided to color the populations S1, S2 with light blue and dark blue, S3, S4 with light green and dark green, S5, S6 with light red and dark red, and S7 with purple. We picked a sequence of colors using the `get_seq_colors()` function. This function returns a sequence of HEX color codes that can be used to specify the coloring scheme for `plot_tree()` function.

```R
colors <- c(get_seq_colors("Blues", 2), get_seq_colors("Greens", 2), get_seq_colors("Reds", 2), get_seq_colors("Purples", 1))
names(colors) <- paste0("S", 1:7)
plot_tree(tree, colors = colors, font_size = 15)
```

![HGDP colorred tree](https://github.com/StoreyLab/superadmixture/blob/main/vignettes/figures/tree-hgdp.png)

We can visualize the estimated admixture proportions using the `barplot_admix` function, which generates the following plot.

```R
barplot_admix(Q_hat, colors = colors, subpops = fam_hgdp$subpop, indiv_on_cols = TRUE)
```

![HGDP admixture proportions](https://github.com/StoreyLab/superadmixture/blob/main/vignettes/figures/admix-props-hgdp.png)

We also can visualize the coancestry among antecedent populations by heatmaps.

```R
par(xpd = TRUE)
heatmap_coanc_antepops(coanc_pops_sup, tl.offset = 1)
```

![HGDP heatmap of coancestry](https://github.com/StoreyLab/superadmixture/blob/main/vignettes/figures/coanc-antepops-hgdp.png)

### Calculating individual-level coancestry under super admixture and standard admixture

We can obtain the corresponding individual-level coancestry under the super admixture model and under the standard admixture model as follows.

```R
coanc_sup <- t(Q_hat) %*% coanc_pops_sup %*% Q_hat
coanc_std <- t(Q_hat) %*% coanc_pops_std %*% Q_hat
```

### Simulating genotypes from the super-admixture model

```R

# approximate ancestral allele frequencies by average allele frequencies
p_anc <- 0.5 * colMeans(X_amr, na.rm = TRUE)

# simulate genotypes according to the double-admixture method
X_sim_amr <- dbl_admixture(p_anc, coanc_pops_sup, Q_hat, geno_only = TRUE)
```

## Citations

* Danfeng Chen, John D. Storey. 2024. "Coancestry superposed on admixed populations yields measures of relatedness at individual-level resolution." bioRxiv doi: [10.1101/2024.12.29.630632](https://www.biorxiv.org/content/10.1101/2024.12.29.630632v1).
* For users of ALStructure for estimating admixture proportions, please cite: Irineo Cabreros, John D. Storey. 2019. "A likelihood-free estimator of population structure bridging admixture models and principal components analysis." doi: [10.1534/genetics.119.302159](https://doi.org/10.1534/genetics.119.302159).
* For users of structured Hardy–Weinberg Equilibrium test, please cite: Wei Hao, John D. Storey. 2019. "Extending tests of Hardy–Weinberg Equilibrium to structured populations." doi: [10.1534/genetics.119.302370](https://doi.org/10.1534/genetics.119.302370). 
