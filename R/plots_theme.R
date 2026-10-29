# One theme for every plot the package draws (diagnostics, fits, comparisons and the report's
# own plots). Monolix's plots in the report keep their own style.

#' Theme of the package's plots
#'
#' `theme_minimal()` without minor gridlines, with bold facet labels on a light
#' grey strip. Groups are coloured with `scale_color_brewer(palette = "Dark2")`.
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme.
#' @export
theme_report <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey95", color = NA))
}

# ggplot() with the package theme; a plot's own theme() tweaks are added after it
ggplot_report <- function(...) ggplot(...) + theme_report()
