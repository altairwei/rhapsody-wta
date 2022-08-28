plotSlingshotCurveOnReduc <- function(sce, reduction = "UMAP", ...) {
  if (!is.null(colnames(reducedDim(sce, reduction))))
    colnames(reducedDim(sce, reduction)) <- NULL

  pseudo.paths <- slingshot::slingPseudotime(sce)
  
  # Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
  # in segments that are shared across paths have similar pseudo-time values in 
  # all paths anyway, so taking the rowMeans is not particularly controversial.
  shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)
  
  # Need to loop over the paths and add each one separately.
  embedded <- slingshot::embedCurves(sce, reduction)
  embedded <- slingshot::slingCurves(embedded)
  curve_data <- lapply(embedded, function(path) data.frame(path$s[path$ord,]))
  
  if (is.data.frame(curve_data))
    curve_list <- list(curve_data)
  else if (is.list(curve_data) && all(sapply(curve_data, is.data.frame)))
    curve_list <- curve_data
  else
    curve_list <- list()
  
  p <- scater::plotReducedDim(
    sce, reduction, colour_by = I(shared.pseudo), ...) +
    ggplot2::coord_fixed()

  for (curve in curve_list)
    p <- p + ggplot2::geom_path(
      data = curve,
      mapping = ggplot2::aes(x = Dim.1, y = Dim.2),
      size = 1.2)
  
  p
}

quickDynSlingshot <- function(sce, reduction = "PCA") {
  dataset <- dynwrap::wrap_expression(
    expression = t(assay(sce, "logcounts")),
    counts = t(assay(sce, "counts"))
  )
  
  dataset <- dynwrap::add_grouping(
    dataset,
    sce$cluster_id
  )
  
  dataset <- dynwrap::add_prior_information(
    dataset,
    dimred = reducedDim(sce, reduction)
  )
  
  model <- dynwrap::infer_trajectory(
    dataset,
    dynmethods::ti_slingshot(),
    give_priors = c("dimred"),
    verbose = TRUE
  )
  
  model
}

runSlingshotAndVelocity <- function(sce, reduction = "PCA", ...) {
  sce <- slingshot::slingshot(
    sce,
    clusterLabels = sce$cluster_id,
    reducedDim = reduction, ...)

  altExp(sce, "scvelo") <- velociraptor::scvelo(
    sce,
    mode = "stochastic",
    use.theirs = TRUE,
    use.dimred = reduction,
    scvelo.params = list(
      filter_and_normalize = list(
        min_shared_counts = 20L, n_top_genes = 2000L),
      moments = list(n_pcs = 30L, n_neighbors = 30L)
    ))

  sce
}

plotVelocityArrow <- function(sce, reduction, ...) {
  velo.out <- altExp(sce, "scvelo")
  sce$velocity_pseudotime <- velo.out$velocity_pseudotime
  embedded <- velociraptor::embedVelocity(reducedDim(sce, reduction), velo.out)
  grid.df <- velociraptor::gridVectors(reducedDim(sce, reduction), embedded)

  scater::plotReducedDim(
      sce, colour_by="velocity_pseudotime",
      dimred = reduction, ...) +
    ggplot2::geom_segment(
      data = grid.df,
      mapping = ggplot2::aes_string(
        x = "start.1",
        y = "start.2", 
        xend = "end.1",
        yend = "end.2"),
      arrow = grid::arrow(
        length = grid::unit(0.02, "inches"),
        type = "closed"),
      size = 0.2) +
    ggplot2::coord_fixed()
}

#' @param sce Slingshot should have been executed.
printDiffTrajectoryReport <- function(sce, heading) {
  # Treatment imbalance
  scores <- condiments::imbalance_score(
    Object = reducedDim(sce, "TSNE"),
    conditions = sce$treatment)
  
  sce$imbalance_score <- scores$scores
  sce$imbalance_scaled_scores <- scores$scaled_scores
  
  df_to_plot <- SingleCellExperiment::reducedDim(sce, "TSNE") |>
    as.data.frame()
  colnames(df_to_plot) <- c("TSNE_1", "TSNE_2")
  df_to_plot$scores <- sce$imbalance_scaled_scores
  pi <- ggplot2::ggplot(df_to_plot,
          ggplot2::aes(x = TSNE_1, y = TSNE_2, col = scores)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(option = "C") +
    ggplot2::coord_fixed() +
    cowplot::theme_cowplot() +
    ggplot2::ggtitle("imbalance score")
  ps <- scater::plotReducedDim(sce, dimred = "TSNE", colour_by = "group_id") +
    ggplot2::coord_fixed()

  old_opts <- knitr::opts_chunk$get()
  knitr::opts_chunk$set(fig.height=7, fig.width=7)
  knitr::knit_print(pi + ps)
  knitr::opts_chunk$restore(old_opts)

  # Slingshot and velocity
  psts <- slingshot::slingPseudotime(sce$slingshot) %>%
    as.data.frame() %>%
    dplyr::mutate(
      cells = rownames(.),
      conditions = sce$treatment) %>%
    tidyr::pivot_longer(dplyr::starts_with("Lineage"),
                        values_to = "pseudotime", names_to = "lineages")

  for (reduc in c("TSNE", "UMAP", "DiffusionMap", "PHATE")) {
    colnames(reducedDim(sce, reduc)) <- NULL
    p1 <- plotSlingshotCurveOnReduc(sce, reduction = reduc) +
      ggplot2::coord_fixed()
    p2 <- plotVelocityArrow(sce, reduction = reduc) +
      ggplot2::coord_fixed()

    old_opts <- knitr::opts_chunk$get()
    knitr::opts_chunk$set(fig.height=6, fig.width=12)
    knitr::knit_print(p1 + p2)
    knitr::opts_chunk$restore(old_opts)
  }

  # Pseudotime distribution
  cat("Kolmogorov-Smirnov Test:")

  MOCK <- dplyr::filter(psts, conditions == "MOCK") |> dplyr::pull(pseudotime)
  PNR2 <- dplyr::filter(psts, conditions == "PNR2") |> dplyr::pull(pseudotime)
  TR4 <- dplyr::filter(psts, conditions == "TR4") |> dplyr::pull(pseudotime)

  cat("\n\n```\n")
  res <- suppressWarnings(ks.test(MOCK, PNR2))
  print(res)
  cat("\n```\n\n")

  cat("\n\n```\n")
  res <- suppressWarnings(ks.test(MOCK, TR4))
  print(res)
  cat("\n```\n\n")

  prog_res <- condiments::progressionTest(
    sce$slingshot, conditions = sce$treatment,
    global = TRUE, lineages = TRUE)
  knitr::kable(prog_res)

  velo.out <- altExp(sce, "scvelo")
  sce$velocity_pseudotime <- velo.out$velocity_pseudotime

  p1 <- ggplot2::ggplot(psts, ggplot2::aes(x = pseudotime, fill = conditions)) +
    ggplot2::geom_density(alpha = .5) +
    ggplot2::scale_fill_brewer(type = "qual") +
    ggplot2::facet_wrap(~lineages) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white")
    )

  p2 <- colData(sce)[, c("velocity_pseudotime", "treatment")] |>
    as.data.frame() |>
    ggplot2::ggplot(ggplot2::aes(x = velocity_pseudotime, fill = treatment)) +
    ggplot2::geom_density(alpha = .5) +
    ggplot2::scale_fill_brewer(type = "qual") +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white")
    )

  old_opts <- knitr::opts_chunk$get()
  knitr::opts_chunk$set(fig.height=6, fig.width=12)
  knitr::knit_print(p1 / p2)
  knitr::opts_chunk$restore(old_opts)
}

blue2white2red <- function(n) {
  pal <- colorRampPalette(c(
    '#081d58', '#081d58', '#253494', '#41b6c4', '#f7fcf0', # Cold
    '#FFFFFF', # White
    '#ffffb2', '#fecc5c', '#e31a1c', '#800026', '#800026' # Warm
  ))
  pal(n)
}

plotPseudotimeHeatmap <- function(mtx,
  palette = colorRamps::matlab.like2, ...) {

  mtx_uni <- unique(mtx)
  if (length(mtx_uni) < 100)
    q <- max(abs(mtx))
  else
    q <- stats::quantile(abs(mtx), probs = .99)
  
  bks <- seq(-q, q, length.out = 9)

  col_fun <- circlize::colorRamp2(
    breaks = bks,
    colors = palette(length(bks))
  )
  
  ComplexHeatmap::Heatmap(
    mtx,
    col = col_fun,
    cluster_columns = FALSE,
    use_raster = TRUE,
    heatmap_legend_param = list(
      title = NULL,
      legend_direction = "vertical",
      legend_height = grid::unit(6, "cm")
    ),
    ...
  )
}

calculateDiffusionMap <- function(
  x, ...,
  exprs_values="logcounts", dimred=NULL, n_dimred=NULL,
  ncomponents = 2, ntop = 500, subset_row = NULL, scale=FALSE,
  transposed=!is.null(dimred)
) {
  mat <- scater:::.get_mat_from_sce(x, exprs_values=exprs_values, dimred=dimred, n_dimred=n_dimred)

  if (!transposed) {
    mat <- scater:::.get_mat_for_reddim(mat, subset_row=subset_row, ntop=ntop, scale=scale) 
  }

  mat <- as.matrix(mat)
  difmap_out <- destiny::DiffusionMap(mat, ...)
  difmap_out@eigenvectors[, seq_len(ncomponents), drop = FALSE]
}

runDiffusionMap <- function(x, ..., altexp=NULL, name="DiffusionMap") {
  if (!is.null(altexp)) {
    y <- altExp(x, altexp)
  } else {
    y <- x
  }
  reducedDim(x, name) <- calculateDiffusionMap(y, ...)
  x 
}

calculatePHATE <- function(
  x, ...,
  exprs_values="logcounts", dimred=NULL, n_dimred=NULL,
  ncomponents = 2, ntop = 500, subset_row = NULL, scale=FALSE,
  transposed=!is.null(dimred)
) {
  # Cell embeddings can be used as input of PHATE, see
  # https://github.com/KrishnaswamyLab/phateR/issues/60
  mat <- scater:::.get_mat_from_sce(x, exprs_values=exprs_values, dimred=dimred, n_dimred=n_dimred)

  if (!transposed) {
    mat <- scater:::.get_mat_for_reddim(mat, subset_row=subset_row, ntop=ntop, scale=scale) 
  }

  # n_cells x n_features
  mat <- as.matrix(mat)

  out <- phateR::phate(mat, ndim = ncomponents, verbose = FALSE, ...)

  out$embedding
}

runPHATE <- function(x, ..., altexp=NULL, name="PHATE") {
  if (!is.null(altexp)) {
    y <- altExp(x, altexp)
  } else {
    y <- x
  }
  reducedDim(x, name) <- calculatePHATE(y, ...)
  x 
}
