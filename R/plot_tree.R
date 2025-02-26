#' Return the node index of the current tree node
#' @param query A string label of the current tree node.
#' @param tr A phylo object from the `ape` package.
#' @return A scalar that corresponds to node index of the current tree node.
#' @export
#' @keywords internal
get_current_node <- function(query, tr) {
  trdata <- tidytree::as_tibble(tr)
  if (query %in% trdata$node) return(query)
  query <- subset(trdata, label==query)$node
  return(query)
}


#' Return the node index of parent tree node
#' @inheritParams get_current_node
#' @return A scalar that corresponds to node index of the parent tree node.
#' @export
#' @keywords internal
get_parent_node <- function(query, tr) {
    query <- get_current_node(query, tr)
    trdata <- tidytree::as_tibble(tr)
    return(tidytree::parent(trdata, query)$node)
}


#' Visualize the coancestry among antecedent populations as a phylogenetic tree
#'
#' @description `plot_tree()` serves as a wrapper function for the \code{\link[ggtree]{ggtree}} function from the ggtree package. The user can use this function to visualize the estimated coancestry matrix \eqn{\Lambda} as a phylogenetic tree.
#'
#' @inheritParams get_current_node
#' @param tip_label An optional string of prefix of the tip label (default to "S").
#' @param font_size An optional float number of font size in pt (default to 14).
#' @param colors An optional named vector of Hex color codes. The names of the vector should have a one-to-one correspondence to the tip labels of `tr`. The values of the vector encoding the colors will be used to plot the tree labels. If no `colors` is provided, this function will return an uncolorred tree.
#'
#' @return A ggplot object of the phylogenetic tree.
#'
#' @note
#' This function depends on the Bioconductor packages `treeio` and `ggtree`. To run this function, the user should install `treeio` and `ggtree` packages mannually.
#' ```
#' if (!require("BiocManager", quietly = TRUE))
#'       install.packages("BiocManager")
#'  BiocManager::install("ggtree")
#'  BiocManager::install("treeio")
#' ```
#'
#' @examples
#' ## Visualizing the coancestry of HGDP data ----------------------------------------
#' data("coanc_pops_hgdp", package = "superadmixture")
#'
#' ## order in the ascending order of coancestry
#' index <- order(diag(coanc_pops_hgdp))
#' coanc_pops_hgdp <- coanc_pops_hgdp[index, index]
#'
#' # label antecedent populations
#' colnames(coanc_pops_hgdp) <- rownames(coanc_pops_hgdp) <- paste0("S", 1:7)
#'
#' ## fit tree
#' tree <- bnpsd::fit_tree(coanc_pops_hgdp)
#'
#' ## plot an uncolorred tree
#' plot_tree(tree, font_size = 15)
#'
#' ## Based on the topology of the tree presented in the previous code chunk,
#' ## we decide to color the populations S1, S2 by light blue and dark blue,
#' ## S3, S4 as light green and dark green, S5, S6 as light red and dark red,
#' ## and the rest with purple (see vignette). We can pick a sequence of colors
#' ## by using the `get_seq_colors()` function and its returned value can be used
#' ## to specify the coloring scheme for `plot_tree()` function.
#' colors <- c(get_seq_colors("Blues", 2),
#'             get_seq_colors("Greens", 2),
#'             get_seq_colors("Reds", 2),
#'             get_seq_colors("Purples", 1))
#' names(colors) <- paste0("S", c(1:4, 5, 6, 7))
#' plot_tree(tree, colors = colors, font_size = 15)
#' @export
#'
plot_tree <- function(tr, tip_label = "S", font_size = 14, colors = NULL) {

    k_antepops <- length(tr$tip.label)

    stopifnot("`tr$tip.label` should be equal to `paste0(tip_label, 1:k_antepops)`" = all(sort(tr$tip.label) == sort(paste0(tip_label, 1:k_antepops))))

    stopifnot("`names(colors)` should match labels of the tree" = is.null(colors) || all(sort(names(colors)) == sort(tr$tip.label)))

    # convert a tr to a tbl
    d  <- treeio::as_tibble(tr) %>%
          dplyr::filter(! is.na(label)) %>%
          dplyr::mutate(subpop = paste0("<i style='font-size:", font_size, "pt'>**", tip_label, "<sub>", gsub(tip_label, "", label), "</sub>**</i>"))

    # plot
    p <-ggtree::"%<+%"(ggtree::ggtree(tr), d) + ggtree::layout_dendrogram() + ggplot2::theme(plot.margin = ggplot2::unit(c(5, 5, 15, 5), "pt"))

    # match colors
    subpop <- dplyr::filter(treeio::as_tibble(tr), !is.na(label))$label

    # color labels if ! is.null(colors)
    if (! is.null(colors)) {
        colors <- colors[match(subpop, names(colors))]
        p <- p + ggtext::geom_richtext(data = ggtree::td_filter(isTip), ggplot2::aes(label = subpop), fill = colors, label.color = NA, alpha = 0.8, hjust = 0.5, vjust = 1)
    } else {
        p <- p + ggtext::geom_richtext(data = ggtree::td_filter(isTip), ggplot2::aes(label = subpop), hjust = 0.5, vjust = 1)
    }
    return(p)
}
