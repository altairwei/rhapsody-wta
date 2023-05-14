#' Preview GGPlot2 Outputs
ggpreview <- function(ggobj, ext = ".png",...) {
  temp.file <- paste0(tempfile(), ext)
  ggplot2::ggsave(filename = temp.file, plot = ggobj, ...)
  utils::browseURL(temp.file)
}

#' Preview Outputs of Python matplotlib
pltpreview <- function(
    fig, ext = ".png",
    width = 7, height = 7, dpi = 300, ...) {
  temp.file <- paste0(tempfile(), ext)
  fig$set_size_inches(width, height)
  fig$savefig(temp.file, dpi = dpi, ...)
  utils::browseURL(temp.file)
}

#' Save Plots Within Powerpoint Slides
ggpowerpoint <- function(ggobj, template, file = NULL, ...) {
  target.file <- file
  if (is.null(file))
    target.file <- paste0(tempfile(), ".pptx")
  
  if (!is.null(template))
    file.copy(template, target.file, overwrite = TRUE)
  
  # TODO: allow create pptx with given size using internal function sof {officer}.
  # Refer to: https://learn.microsoft.com/en-us/office/open-xml/open-xml-sdk
  margins = c(top = 0.5, right = 0.5,bottom = 0.5, left = 0.5)
  doc <- officer::read_pptx(path = target.file)
  doc = officer::add_slide(doc, layout = "Blank", master = "Office Theme")
  pagesize = export:::get.slide.size(doc)
  pagesize["width"] = pagesize["width"] - (margins["left"] + margins["right"])
  pagesize["height"] = pagesize["height"] - (margins["top"] + margins["bottom"])
  
  doc = officer::ph_with(
    x = doc,
    value = rvg::dml(ggobj = ggobj, ...),
    location = officer::ph_location(
      left = margins["left"], top = margins["top"],
      width = pagesize["width"], height = pagesize["height"]))
  
  print(doc, target = target.file)
  
  if (is.null(file))
    utils::browseURL(target.file)
}

# Preview Grid Plots, such ComlexHeatmap
plotpreview <- function(x, width = 7, height = 7, res = 300, ...) {
  temp.file <- paste0(tempfile(), ".png")
  png(filename = temp.file, units = "in",
      width = width, height = height, res = res, ...)
  print(x)
  dev.off()
  utils::browseURL(temp.file)
  invisible()
}

rotate_x_labels <- function(angle = 45, hjust = 1, ...) {
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = angle, hjust = hjust, ...))
}

verticalize_x_labels <- function(angle = 90, hjust = 1, vjust = 0.5, ...) {
  rotate_x_labels(angle = angle, hjust = hjust, vjust = vjust, ...)
}

center_plot_title <- function(...) {
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, ...))
}

remove_axis <- function(...) {
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    ...)
}

remove_x_axis <- function(...) {
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_blank(),
    ...)
}

remove_legend <- function(...) {
  ggplot2::theme(
    legend.position = "none", validate = TRUE, ...)
}

remove_legend_title <- function(...) {
  ggplot2::theme(
    legend.title = ggplot2::element_blank(), ...)
}

empty_strip <- function(...) {
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    ...
  )
}

remove_strip <- function(...) {
  ggplot2::theme(
    strip.text = ggplot2::element_blank() , 
    strip.background = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0) , units = "lines"),
    ...)
}

legend_override <- function(key, values, ...) {
  args <- list()
  args[[key]] <- ggplot2::guide_legend(
    override.aes = values, ...)
  do.call(ggplot2::guides, args)
}

font_plot_tag <- function(size = 11, ...) {
  ggplot2::theme(plot.tag = ggplot2::element_text(
    size = size, face = "bold", ...))
}

font_text_title <- function(...) {
  ggplot2::theme(text = ggplot2::element_text(...),
                 title = ggplot2::element_text(...))
}

theme_dimred <- function(
    xlength = 0.3, ylength = 0.3, 
    arrow = grid::arrow(
      angle = 15, length = grid::unit(0.15, "inches"),
      type = "closed"),
    ...) {
  ggplot2::theme(
    axis.line = ggplot2::element_blank(), 
    axis.ticks = ggplot2::element_blank(), 
    axis.text = ggplot2::element_blank(),
    axis.line.x.bottom = tidydr::element_line2(
      id = 1, xlength = xlength, arrow = arrow),
    axis.line.y.left = tidydr::element_line2(
      id = 2, ylength = ylength, arrow = arrow), 
    axis.title = ggplot2::element_text(hjust=0.1),
    ...
  )
}

theme_upset <- function(...) {
  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(
      colour = "black", linewidth = ggplot2::rel(1)),
    legend.key = ggplot2::element_blank(), 
    strip.background = ggplot2::element_rect(
      fill = "white", colour = "black", linewidth = ggplot2::rel(2)),
    ...
  )
}

#' Generate Radar Plot From GSVA Data Frame 
plotRadarGSVA <- function(gsva, pathws) {
  dplyr::filter(gsva,
                Pathway %in% pathws) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = cellType,
        y = ES,
        group = treatment,
        color = treatment,
        fill = treatment
      )) +
    ggplot2::geom_polygon(linewidth = 1, alpha = .1) +
    ggplot2::geom_point() +
    ggplot2::scale_color_brewer(type = "qual", palette = "Dark2") +
    ggplot2::scale_fill_brewer(type = "qual", palette = "Dark2") +
    see::coord_radar() +
    ggplot2::facet_wrap(~ time, nrow = 1) +
    see::theme_radar() +
    empty_strip() +
    NULL
}

#' Generate Line Plot with Facets From GSVA Data Frame
plotLineFacetGSVA <- function(gsva, pathws) {
  df <- dplyr::filter(gsva, Pathway %in% pathws)
  
  df_0 <- dplyr::filter(df, time == "0DPI")
  
  df |>
    dplyr::add_row(dplyr::mutate(df_0, treatment = "PNR2")) |>
    dplyr::add_row(dplyr::mutate(df_0, treatment = "TR4")) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = time,
        y = ES,
        group = treatment,
        color = treatment,
        fill = treatment
      )) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point() +
    ggplot2::scale_color_brewer(type = "qual", palette = "Dark2") +
    ggplot2::facet_wrap(~ cellType, nrow = 1) +
    ggplot2::theme_bw() +
    NULL
}

plotDEGsHeatmap <- function(
    mtx, genes, colAnno,
    max = 5, asterisk = TRUE,
    show_row_names = TRUE, draw = TRUE, ...) {
  
  detected <- intersect(rownames(mtx), genes)
  mtx <- mtx[detected,]
  
  if (length(detected) != length(genes))
    warning(sprintf(
      "Only %i out of %i genes were detected in DEGs\n",
      length(detected), length(genes)))
  
  outliers_idx <- which(abs(mtx) > max, arr.ind = TRUE)
  if (length(outliers_idx) != 0)
    warning(sprintf(paste0("%i value%s with abs(logFC) greater than %i were shown",
                           " in maximum color and marked with white asterisk.\n"),
                    nrow(outliers_idx),
                    if (nrow(outliers_idx) > 1) "s" else "",
                    max))
  
  p <- ComplexHeatmap::Heatmap(
    matrix = mtx,
    name = "Log2 Fold Change",
    col = circlize::colorRamp2(c(-max, 0, max), c("blue", "white", "red")),
    use_raster = TRUE,
    na_col = "grey",
    cluster_rows = function(x) {
      dend_mtx <- x
      dend_mtx[is.na(dend_mtx)] <- 0
      dendsort::dendsort(
        fastcluster::hclust(dist(dend_mtx)),
        isReverse = TRUE)
    },
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = show_row_names,
    top_annotation = colAnno,
    heatmap_legend_param = list(
      title_position = "leftcenter-rot",
      legend_height  = grid::unit(4, "cm")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- abs(mtx[i, j])
      if (asterisk && !is.na(val) && val > max)
        grid::grid.text("★", x, y, gp = grid::gpar(
          fontsize = 8, col = "white", fontface = "bold"))
    },
    ...)

  if (draw)
    ComplexHeatmap::draw(p, merge_legend = TRUE)
  else
    p
}
