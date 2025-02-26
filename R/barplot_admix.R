#' Visualize admixture proportions of each individual
#'
#' @description
#' `barplot_admix()` serves as a wrapper function for the
#' \code{\link[ggplot2]{geom_bar}} function from the ggplot2 package
#' to visualize the admixture proportion.
#'
#' @param admix_props A matrix of admixture proportions. The dimension of `admix_props` can be \eqn{K} (number of antecedent populations) \eqn{\times} \eqn{n} (number of individuals) or \eqn{n \times K}. Set `indiv_on_cols` to be `FALSE` if the dimension of `admix_props` is \eqn{n \times K}.
#' @param colors An optional vector of Hex color codes. The length of `colors` should be equal to \eqn{K}. Each hex color will be used to plot the corresponding antecedent populations. If no `colors` is provided, this function will use the Accent palette as default.
#' @param subpops An optional vector of strings each representing population label of the corresponding individual. If `subpops` is set, `admix_props` should be arranged in a way that individuals from the same population are grouped together. For example, `admix_props` has 7 individuals, two from "S1", two from "S2", the rest from "S3", then `subpops` can be `c("S1", "S1", "S2", "S2", "S3", "S3", "S3")`.
#' @param indiv_on_cols An optional boolean that specifies the dimension of the input `admix_props`. If `FALSE`, `admix_props` has individuals on rows and populations on cols; if `TRUE` (default), populations are on rows and individuals on columns.
#' @param base_size An optional integer of base font size for \code{\link[ggplot2]{geom_bar}} (default = 14).
#'
#' @return A ggplot object
#'
#' @examples
#' ## Visualizing admixture proportions of HGDP data ----------------------------------------
#'
#' ## load admixture proportions, population-level coancestry and FAM data
#' data("coanc_pops_hgdp",  package = "superadmixture")
#' data("admix_props_hgdp", package = "superadmixture")
#' data("fam_hgdp", package = "superadmixture")
#'
#' ## order in the ascending order of coancestry
#' index <- order(diag(coanc_pops_hgdp))
#' admix_props_hgdp <- admix_props_hgdp[index, ]
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
#' names(colors) <- paste0("S", 1:7)
#'
#' ## plot admixture proportions without population labels
#' barplot_admix(admix_props_hgdp, colors = colors, indiv_on_cols = TRUE)
#'
#' ## plot admixture proportions with population labels
#' barplot_admix(admix_props_hgdp, colors = colors, subpops = fam_hgdp$subpop, indiv_on_cols = TRUE)
#' @export
barplot_admix <- function(admix_props, colors = NULL, subpops = NULL, indiv_on_cols = TRUE, base_size = 14) {

    # transpose if needed
    check_admix_proportions(admix_props)
    if (indiv_on_cols) admix_props <- t(admix_props)

    # check whether `length(colors) == ncol(admix_props)`
    stopifnot("`length(colors)` should be equal to `ncol(admix_props)`" = is.null(colors) || length(colors) == ncol(admix_props))

    # get population labels
    k_subpops <- ncol(admix_props)
    subpop_label <- paste0("S", 1:k_subpops)

    # convert the admixture proportions to the long format
    admix_props_long <- tibble::as_tibble(admix_props) %>%
        dplyr::rename_at(dplyr::vars(paste0("V", 1:k_subpops)), ~subpop_label) %>%
        dplyr::mutate(idx = 1:dplyr::n()) %>%
        tidyr::pivot_longer(!idx, names_to="subpop", values_to="admix_prop") %>%
        dplyr::mutate(subpop=factor(subpop, subpop_label))

    # create names for colors
    if (!is.null(names(colors))) {
        if (any(sort(as.character(names(colors))) != sort(as.character(subpop_label))))
            stop("`names(colors)` should match `subpop_label`")
    } else if (! is.null(colors)) {
        names(colors) <- subpop_label
    } else {
        colors <- get_seq_colors("Accent", k_subpops)
        names(colors) <- subpop_label
    }

    # create labels
    labels <- gsub("S", "S[", names(colors))
    labels <- paste0(labels, "]")
    names(labels) <- names(colors)

    # re-order colors and labels
    idx <- match(subpop_label, names(colors))
    colors <- colors[idx]
    labels <- labels[idx]

    # plot the admixture props
    p <- ggplot2::ggplot(admix_props_long, ggplot2::aes_string(x="idx", y="admix_prop", col="subpop", fill="subpop")) +
        ggplot2::geom_bar(stat='identity') +
        ggplot2::scale_fill_manual( labels = parse(text = labels), values = colors) +
        ggplot2::scale_color_manual(labels = parse(text = labels), values = colors, guide = "none") +
        ggplot2::scale_x_continuous("Individuals", expand=c(0, 0)) +
        ggplot2::scale_y_continuous("Admix Props", expand=c(0, 0)) +
        ggplot2::theme_bw(base_size=base_size) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            legend.title = ggplot2::element_blank(),
            axis.text    = ggplot2::element_blank(),
            axis.ticks   = ggplot2::element_blank()) +
        ggplot2::guides(fill =ggplot2::guide_legend(title="Sub-pops"))


    if (! is.null(subpops)) {
        if (! is.factor(subpops)) subpops <- factor(subpops, unique(subpops))
        xintercept <- cumsum(table(subpops))
        breaks <- xintercept - 0.5 * as.numeric(table(subpops))
        p <- p + ggplot2::scale_x_continuous(
                "Individuals",
                expand=c(0, 0),
                breaks=breaks,
                labels=unique(subpops)) +
            ggplot2::geom_vline(
                xintercept=utils::head(xintercept,-1),
                linetype="solid",
                linewidth=1) +
            ggplot2::theme(
                axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    }
    # return
    return(p)

}
