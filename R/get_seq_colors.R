#' Get a sequence of colors from RColorBrewer palettes
#'
#' @description
#' `get_seq_colors()` returns a sequence of Hex color codes from RColorBrewer package.
#'
#' @param palette A string of RColorBrewer palette name.
#' @param n_colors An integer of number of colors. `n_colors` must be greater than 0.
#'
#' @return A vector of Hex color codes.
#'
#' @note
#' `palette` must be one of the following: BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral, Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3, Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd.
#'
#' @examples
#' ## Getting three colors from 'Reds' palette --------------------------
#' get_seq_colors("Reds", 3)
#'
#' ## Getting 20 colors from 'Spectral' palette -------------------------
#' get_seq_colors("Spectral", 20)
#'
#' @export
#' @keywords internal
get_seq_colors <- function(palette, n_colors) {
    stopifnot("`n_colors` must be greater than 0" = n_colors > 0)
    if (! palette %in% rownames(RColorBrewer::brewer.pal.info)){
      stop(paste0("Invalid palette name! Please specify a palette from the following: ",
                  paste(rownames(RColorBrewer::brewer.pal.info), collapse = ", "), "."))
    }
    if (n_colors < 3) {
        return(RColorBrewer::brewer.pal(4, name = palette)[(4 - n_colors): 3])
    } else if (n_colors < 4) {
        return(RColorBrewer::brewer.pal(4, name = palette)[(5 - n_colors): 4])
    } else if (n_colors > 9) {
        return(grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, name = palette))(n_colors))
    } else {
        return(grDevices::colorRampPalette(RColorBrewer::brewer.pal(n_colors, name = palette))(n_colors))
    }
}
