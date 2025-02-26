#' Visualize the coancestry matrix between antecedent populations as a heatmap
#'
#' @description `heatmap_coanc_antepops()` serves as a wrapper function for the \code{\link[corrplot]{corrplot}} function from the corrplot package to plot the coancestry between antecedent population \eqn{\boldsymbol{\Lambda}} as a heatmap.
#'
#' @inheritParams est_paras_dbladmix
#' @inheritParams plot_tree
#' @param tl.cex Numeric, for the size of text label (variable names).
#' @param cl.cex Numeric, cex of number-label in color-legend.
#' @param tl.offset Numeric, for text label, see \code{\link[corrplot]{corrplot}}.
#'
#' @return A heatmap of \eqn{\boldsymbol{\Lambda}}.
#'
#' @examples
#' ## Visualizing the coancestry of HGDP data ----------------------------------------
#' data("coanc_pops_hgdp", package = "superadmixture")
#'
#' ## order in the ascending order of coancestry
#' index <- order(diag(coanc_pops_hgdp))
#' coanc_pops_hgdp <- coanc_pops_hgdp[index, index]
#'
#' ## label antecedent populations
#' colnames(coanc_pops_hgdp) <- rownames(coanc_pops_hgdp) <- paste0("S", 1:7)
#'
#' ## plot the heatmap of coancestry
#' par(xpd = TRUE)
#' heatmap_coanc_antepops(coanc_pops_hgdp)
#'
#' @export
heatmap_coanc_antepops <- function(coanc_antepops, tip_label = "S", tl.cex = 1.3, cl.cex = 1, tl.offset = 0.6) {

  k_antepops <- nrow(coanc_antepops)

    stopifnot("colnames(coanc_antepops) should be equal to `paste0(tip_label, 1:k_antepops)`" = all(colnames(coanc_antepops) == paste0(tip_label, 1:k_antepops)))
    stopifnot("rownames(coanc_antepops) should be equal to `paste0(tip_label, 1:k_antepops)`" = all(rownames(coanc_antepops) == paste0(tip_label, 1:k_antepops)))

    colnames(coanc_antepops) <- gsub(tip_label, paste0(":", tip_label, "["), colnames(coanc_antepops))
    colnames(coanc_antepops) <- paste0(colnames(coanc_antepops), "]")
    rownames(coanc_antepops) <- gsub(tip_label, paste0(":", tip_label, "["), rownames(coanc_antepops))
    rownames(coanc_antepops) <- paste0(rownames(coanc_antepops), "]")

    col <- grDevices::colorRampPalette(c("#FFFFFF", "#EE9988", "#BB4444"))
    corrplot::corrplot(coanc_antepops, method = "color", col=col(200),
         type = "lower", col.lim = c(0, max(coanc_antepops) + 0.05), is.cor = FALSE,
         font = 2, tl.cex = tl.cex, tl.col = "black", tl.srt = 0, tl.offset = tl.offset,
         cl.cex = cl.cex, cl.length = 5, addgrid.col = "white", diag = TRUE,  mar = c(0.1, 1, 1, 0.1))

}
