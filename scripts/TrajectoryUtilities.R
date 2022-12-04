suppressMessages(library(SingleCellExperiment))
suppressMessages(library(magrittr))

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

runUpstream <- function(sce, ...) {
  hvgs <- scran::getTopHVGs(sce, n = 2000)
  sce <- scran::fixedPCA(sce, subset.row = hvgs)
  sce$sample <- as.character(sce$sample_id)
  sce <- harmony::RunHarmony(
    sce,
    ...,
    kmeans_init_nstart=20, kmeans_init_iter_max=1000,
    epsilon.cluster=-Inf, epsilon.harmony=-Inf,
    plot_convergence = TRUE, verbose = FALSE)
  sce
}

runSlingshotAndVelocity <- function(sce, reduction = "PCA", ...) {
  sce <- slingshot::slingshot(
    sce,
    clusterLabels = sce$cluster_id,
    reducedDim = reduction, ...)

  basilisk::setBasiliskShared(FALSE)
  # FIXME: set `enforce=True` to normalize X which is corrected for
  # ambient RNA, otherwise use `assay.X="logcounts"` instead.
  altExp(sce, "scvelo") <- velociraptor::scvelo(
    sce,
    mode = "stochastic",
    use.theirs = TRUE,
    use.dimred = reduction,
    scvelo.params = list(
      filter_and_normalize = list(
        min_shared_counts = 20L, n_top_genes = 2000L)
    ))

  sce
}

runTradeSeq <- function(sce, subset_row = NULL, nknots = 6, ...) {
  altExp(sce, "tradeSeq") <- tradeSeq::fitGAM(
    # fitGAM needs counts rather than logcounts, see
    # the source code of SCE-version fitGAM.
    counts = sce,
    conditions = factor(sce$treatment),
    nknots = nknots,
    sce = TRUE, verbose = TRUE,
    genes = if (!is.null(subset_row))
      subset_row else seq_len(nrow(sce)),
    ...
  )

  sce
}

runScVelo <- function(x,
  dimred = NULL, use.dimred = NULL,
  assay.X="counts", assay.spliced="spliced", assay.unspliced="unspliced",
  mode = c('steady_state', 'deterministic', 'stochastic', 'dynamical'),
  scvelo.params = list(), ncomponents = 30, return_sce = FALSE) {

  if (is.null(dimred)) {
    if (!is.null(use.dimred)) {
      dimred <- reducedDim(x, use.dimred)
    }
  }

  spliced <- assay(x, assay.spliced)
  unspliced <- assay(x, assay.unspliced)
  X <- assay(x, assay.X)

  refdim <- as.integer(dim(spliced))
  if (!identical(refdim, as.integer(dim(unspliced))) || !identical(refdim, as.integer(dim(X)))) {
    stop("matrices in 'x' must have the same dimensions")
  }

  X <- t(X)
  spliced <- t(spliced)
  unspliced <- t(unspliced)

  and <- reticulate::import("anndata")
  scv <- reticulate::import("scvelo")
  adata <- and$AnnData(X, layers=list(spliced=spliced, unspliced=unspliced))
  adata$obs_names <- rownames(spliced)
  adata$var_names <- colnames(spliced)

  if (!is.null(dimred)) {
    adata$obsm <- list(X_pca = dimred)
  }

  # Preprocessing requisites consist of
  #  1. gene selection by detection (with a minimum number of counts) and high variability (dispersion)
  #  2. normalizing every cell by its total size and logarithmizing X.
  #  3. Filtering and normalization is applied in the same vein to spliced/unspliced counts and X.
  #  4. Logarithmizing is only applied to X.
  #  5. If X is already preprocessed from former analysis, it will not be touched.
  do.call(scv$pp$filter_and_normalize, c(list(data=adata), scvelo.params$filter_and_normalize))

  do.call(scv$pp$moments, c(list(data=adata), scvelo.params$moments))

  if (mode=="dynamical") {
    do.call(scv$tl$recover_dynamics, c(list(data=adata), scvelo.params$recover_dynamics))
  }

  scvelo.params$velocity$mode <- mode
  do.call(scv$tl$velocity, c(list(data=adata), scvelo.params$velocity))

  do.call(scv$tl$velocity_graph, c(list(data=adata), scvelo.params$velocity_graph))

  do.call(scv$tl$velocity_pseudotime, c(list(adata=adata), scvelo.params$velocity_pseudotime))

  if (mode=="dynamical") {
    do.call(scv$tl$latent_time, c(list(data=adata), scvelo.params$latent_time))
  }

  do.call(scv$tl$velocity_confidence, c(list(data=adata), scvelo.params$velocity_confidence))

  if (return_sce)
    return(AnnData2SCE(adata))
  else
    return(adata)
}

runDynamoVelocity <- function(x,
  dimred = NULL, use.dimred = NULL,
  assay.X="counts", assay.spliced="spliced", assay.unspliced="unspliced",
  mode = c("deterministic", "stochastic", "negbin"),
  n_top_genes = 2000L, py.params = list(),
  ...) {

  scv <- reticulate::import("scvelo")
  and <- reticulate::import("anndata")
  dyn <- reticulate::import("dynamo")

  if (is.null(dimred)) {
    if (!is.null(use.dimred)) {
      dimred <- reducedDim(x, use.dimred)
    }
  }

  spliced <- assay(x, assay.spliced)
  unspliced <- assay(x, assay.unspliced)
  X <- assay(x, assay.X)

  refdim <- as.integer(dim(spliced))
  if (!identical(refdim, as.integer(dim(unspliced))) || !identical(refdim, as.integer(dim(X)))) {
    stop("matrices in 'x' must have the same dimensions")
  }

  X <- t(X)
  spliced <- t(spliced)
  unspliced <- t(unspliced)

  and <- reticulate::import("anndata")
  scv <- reticulate::import("scvelo")
  adata <- and$AnnData(X, layers=list(spliced=spliced, unspliced=unspliced))
  adata$obs_names <- rownames(spliced)
  adata$var_names <- colnames(spliced)

  preprocessor <- dyn$preprocessing$Preprocessor(
    convert_gene_name_function = NULL
  )
  preprocessor$config_seurat_recipe(adata)
  preprocessor$standardize_adata(adata = adata, tkey = NULL, experiment_type = "conventional")
  preprocessor$filter_genes_by_outliers(adata)
  preprocessor$normalize_by_cells(adata)
  preprocessor$select_genes_kwargs <- list(
    "recipe" = "seurat", "n_top_genes" = as.integer(n_top_genes))
  preprocessor$select_genes(adata)
  preprocessor$log1p(adata)

  if (!is.null(dimred)) {
    adata$obsm <- list(X_pca = dimred)
  }

  # Experiment type
  # adata$uns$update(list(
  #   "pp" = list(
  #     has_splicing = TRUE,
  #     has_labeling = FALSE,
  #     splicing_labeling = FALSE,
  #     has_protein = FALSE,
  #     experiment_type = "conventional",
  #     experiment_layers = c("X", "spliced", "unspliced"),
  #     experiment_total_layers = NULL,
  #     tkey = NULL,
  #     norm_method = NULL
  #   )
  # ))

  mode <- match.arg(mode)
  do.call(dyn$tl$moments, c(list(adata = adata, X_data = adata$obsm['X_pca']), py.params$moments))
  do.call(dyn$tl$dynamics, c(list(adata = adata, model = mode), py.params$dynamics))

  adata
}

runUniTVelo <- function(x, dimred = NULL, use.dimred = NULL,
  assay.X="counts", assay.spliced="spliced", assay.unspliced="unspliced",
  mode = c("unified-time", "independent"), GPU = -1L, label = "cellType", normalize=TRUE,
  ncomponents = 30, return_sce = FALSE
) {
  if (is.null(dimred)) {
    if (!is.null(use.dimred)) {
      dimred <- reducedDim(x, use.dimred)
    }
  }

  spliced <- assay(x, assay.spliced)
  unspliced <- assay(x, assay.unspliced)
  X <- assay(x, assay.X)

  refdim <- as.integer(dim(spliced))
  if (!identical(refdim, as.integer(dim(unspliced))) || !identical(refdim, as.integer(dim(X)))) {
    stop("matrices in 'x' must have the same dimensions")
  }

  X <- t(X)
  spliced <- t(spliced)
  unspliced <- t(unspliced)

  and <- reticulate::import("anndata")
  adata <- and$AnnData(X, layers=list(spliced=spliced, unspliced=unspliced))
  adata$obs_names <- rownames(spliced)
  adata$var_names <- colnames(spliced)
  adata$obs$cellType <- colData(x)[[label]]

  if (!is.null(dimred)) {
    adata$obsm <- list(X_pca = dimred)
  }

  utv <- reticulate::import("unitvelo")
  mode <- match.arg(mode)
  velo <- utv$config$Configuration()
  velo$R2_ADJUST <- TRUE
  velo$IROOT <- NULL
  velo$FIT_OPTION <- if (mode == "unified-time") '1' else '2'
  velo$GPU <- GPU
  velo$BASIS <- "pca"

  do.call(utv$run_model, list(
    adata=adata, label="cellType",
    config_file=velo, normalize=normalize))

  if (return_sce)
    return(AnnData2SCE(adata))
  else
    return(adata)
}

runScanorama <- function(sce, nfeatures = 2000L, ..., convert = TRUE) {
  obj <- Seurat::CreateSeuratObject(
    counts = assay(sce, "counts"),
    assay = "RNA",
    meta.data = data.frame(
      row.names = colnames(sce),
      sample_id = sce$sample_id,
      cell_type = sce$cluster_id)
  )

  obj_list <- Seurat::SplitObject(obj, split.by = "sample_id")
  datasets <- lapply(obj_list, function(obj) {
    obj %>%
      Seurat::NormalizeData(verbose = FALSE) %>%
      Seurat::FindVariableFeatures(
        selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% 
      Seurat::ScaleData(verbose = FALSE) %>%
      Seurat::GetAssayData(assay = "RNA", slot = "scale.data") %>%
      t()
  })
  names(datasets) <- NULL
  genes <- lapply(datasets, function(data) colnames(data))

  scanorama <- reticulate::import("scanorama")
  results <- scanorama$correct(
    datasets, genes, return_dimred = TRUE,
    dimred = 30L, ds_names = names(obj_list))

  results
}

plotSlingshotCurveOnReduc <- function(sce, dimred = "UMAP", path_size = 3, ...) {
  if (!is.null(colnames(reducedDim(sce, dimred))))
    colnames(reducedDim(sce, dimred)) <- NULL

  pseudo.paths <- slingshot::slingPseudotime(sce)

  # Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
  # in segments that are shared across paths have similar pseudo-time values in 
  # all paths anyway, so taking the rowMeans is not particularly controversial.
  shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

  # Need to loop over the paths and add each one separately.
  dots <- list(...)
  if ("ncomponents" %in% names(dots)) {
    ncomps <- if (length(dots$ncomponents) == 1)
      seq_len(dots$ncomponents) else dots$ncomponents
    dimred <- reducedDim(sce, dimred)[, ncomps]
  }

  embedded <- slingshot::embedCurves(sce, dimred)
  embedded <- slingshot::slingCurves(embedded)
  curve_data <- lapply(embedded, function(path) data.frame(path$s[path$ord,]))

  if (is.data.frame(curve_data))
    curve_list <- list(curve_data)
  else if (is.list(curve_data) && all(sapply(curve_data, is.data.frame)))
    curve_list <- curve_data
  else
    curve_list <- list()


  p <- scater::plotReducedDim(
    sce, dimred = dimred, colour_by = I(shared.pseudo), ...)

  for (curve in curve_list)
    p <- p + ggplot2::geom_path(
      data = curve,
      mapping = ggplot2::aes(x = Dim.1, y = Dim.2),
      size = path_size)

  p
}

plotLineageCurveOnReduc <- function(sce, lineage, dimred = "UMAP", ncomponents = 2, ...) {
  stopifnot(is.numeric(lineage))
  
  if (length(ncomponents) == 1L) {
    to_plot <- seq_len(ncomponents)
  }
  else {
    to_plot <- ncomponents
  }

  embeddings <- reducedDim(sce, dimred)[, to_plot]
  colnames(embeddings) <- paste0(dimred, ".", to_plot)
  embedded <- slingshot::slingCurves(
    slingshot::embedCurves(sce, embeddings))
  pseudo.paths <- slingshot::slingPseudotime(sce)
  
  scater::plotReducedDim(sce,
                         dimred = dimred, ncomponents = ncomponents,
                         colour_by = I(pseudo.paths[, lineage]), ...) +
    ggplot2::geom_path(
      data = data.frame(embedded[[lineage]]$s[embedded[[lineage]]$ord,]),
      mapping = ggplot2::aes_string(
        x = colnames(embeddings)[1],
        y = colnames(embeddings)[2]),
      size = 1.2)
}

plotLineagesOnReduc <- function(
    sce, dimred = "UMAP", ncomponents = 2,
    path_size = 3, ...) {
  if (length(ncomponents) == 1L)
    to_plot <- seq_len(ncomponents)
  else
    to_plot <- ncomponents
  
  embeddings <- reducedDim(sce, dimred)[, to_plot]
  colnames(embeddings) <- paste0(dimred, ".", to_plot)
  
  pseudo.paths <- slingshot::slingPseudotime(sce)
  embedded <- slingshot::slingCurves(
    slingshot::embedCurves(sce, embeddings))
  names(embedded) <- colnames(pseudo.paths)

  piclist <- purrr::imap(embedded, function(emb, name) {
    scater::plotReducedDim(
      sce, dimred = dimred, ncomponents = ncomponents,
      colour_by = I(pseudo.paths[, name]), ...) +
      ggplot2::geom_path(
        data = data.frame(emb$s[emb$ord,]),
        mapping = ggplot2::aes_string(
          x = colnames(embeddings)[1],
          y = colnames(embeddings)[2]),
        size = path_size) +
      ggplot2::ggtitle(name) +
      ggplot2::coord_fixed()
  })
  
  patchwork::wrap_plots(piclist)
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
  cat(heading, "Treatment imbalance\n\n")
  scores <- condiments::imbalance_score(
    Object = reducedDim(sce, "TSNE"),
    conditions = sce$treatment)
  
  sce$imbalance_score <- scores$scores
  sce$imbalance_scaled_scores <- scores$scaled_scores
  
  df_to_plot <- SingleCellExperiment::reducedDim(sce, "TSNE") |>
    as.data.frame()
  colnames(df_to_plot) <- c("TSNE_1", "TSNE_2")
  df_to_plot$scores <- sce$imbalance_scaled_scores
  pi <- suppressWarnings(
      ggplot2::ggplot(df_to_plot,
            ggplot2::aes(x = TSNE_1, y = TSNE_2, col = scores)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_viridis_c(option = "C") +
      ggplot2::coord_fixed() +
      cowplot::theme_cowplot() +
      ggplot2::ggtitle("imbalance score")
  )
  ps <- suppressWarnings(
      scater::plotReducedDim(sce, dimred = "TSNE", colour_by = "group_id") +
        ggplot2::coord_fixed()
  )

  print(pi + ps)

  cat("\n\n")

  cat(heading, "Slingshot and velocity {.tabset}\n\n")
  psts <- slingshot::slingPseudotime(sce$slingshot) %>%
    as.data.frame() %>%
    dplyr::mutate(
      cells = rownames(.),
      conditions = sce$treatment) %>%
    tidyr::pivot_longer(dplyr::starts_with("Lineage"),
                        values_to = "pseudotime", names_to = "lineages")

  for (reduc in c("TSNE", "UMAP", "DiffusionMap", "PHATE")) {
    cat(paste0(heading, "#"), reduc, "\n\n")

    colnames(reducedDim(sce, reduc)) <- NULL
    p1 <- suppressMessages(
      plotSlingshotCurveOnReduc(sce, reduction = reduc) +
      ggplot2::coord_fixed()
    )
    p2 <- suppressMessages(
      plotVelocityArrow(sce, reduction = reduc) +
      ggplot2::coord_fixed()
    )

    print(p1 + p2)

    cat("\n\n")
  }

  cat(heading, "Pseudotime distribution\n\n")

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

  print(p1 / p2)

  cat("\n\n")

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

  prog_res <- suppressMessages(
    condiments::progressionTest(
      sce$slingshot, conditions = sce$treatment,
      global = TRUE, lineages = TRUE)
  )
  print(knitr::kable(prog_res))
}

blue2white2red <- function(n) {
  pal <- colorRampPalette(c(
    '#081d58', '#081d58', '#253494', '#41b6c4', '#f7fcf0', # Cold
    '#FFFFFF', # White
    '#ffffb2', '#fecc5c', '#e31a1c', '#800026', '#800026' # Warm
  ))
  pal(n)
}

plotPseudotimeHeatmap <- function(
    mtx, cluster_rows = TRUE, query_set = FALSE,
    seriation = FALSE, seriation_method = "GW",
    palette = colorRamps::matlab.like2,
    color_branches = FALSE,
    dend_k = NULL,
    ...) {

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

  row_dist <- as.dist(1 - cor(t(mtx)))
  row_dist[is.na(row_dist)] <- 1

  if (seriation) {
    o1 <- seriation::seriate(row_dist, method = seriation_method)
    # Currently, the permutation vector can be stored as a simple
    # integer vector or as an object of class hclust.
    row_dend <- as.dendrogram(o1[[1]])
  } else {
    row_dend <- hclust(row_dist, method = "ward.D2")
    if (color_branches)
      row_dend <- dendextend::color_branches(row_dend, k = dend_k)
  }

  if (cluster_rows) {
    show_row_dend <- TRUE
    cluster_rows <- row_dend
  } else {
    cluster_rows <- FALSE
    show_row_dend = FALSE
  }

  ph_res <- ComplexHeatmap::Heatmap(
    mtx,

    # Remove name from fill legend
    name = NULL,
    use_raster = TRUE,

    # Keep original row/col order
    row_order = rownames(mtx),
    column_order = colnames(mtx),
    col = col_fun,

    # Add left annotation (legend with tumor/normal)
    # right_annotation = gene_ann,
    # ACTUAL SPLIT by sample group
    cluster_rows = cluster_rows,
    row_split = dend_k,
    show_row_names = FALSE,
    show_row_dend = show_row_dend,

    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,

    heatmap_legend_param = list(
      title = NULL,
      legend_direction = "vertical",
      legend_height = grid::unit(6, "cm")
    ),
    ...
  )

  if (query_set) {
    query_results <- list(heatmap = ph_res, row_dend = row_dend)
    return(query_results)
  } else {
    return(ph_res)
  }
}

calculateDiffusionMap <- function(
  x, ...,
  exprs_values="logcounts", dimred=NULL, n_dimred=NULL,
  ncomponents = 2, ntop = 500, subset_row = NULL, scale=FALSE,
  transposed=!is.null(dimred), seed = 1L
) {
  mat <- scater:::.get_mat_from_sce(x, exprs_values=exprs_values, dimred=dimred, n_dimred=n_dimred)

  if (!transposed) {
    mat <- scater:::.get_mat_for_reddim(mat, subset_row=subset_row, ntop=ntop, scale=scale) 
  }

  mat <- as.matrix(mat)
  set.seed(seed)
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

wrapSlingshotToDynverse <- function(sds) {
  cluster <- slingshot::slingClusterLabels(sds)
  start_cell <- apply(slingshot::slingPseudotime(sds), 1, min) |>
    sort() |> head(1) |> names()
  start.clus <- names(which(cluster[start_cell,] == 1))

  # satisfy r cmd check
  from <- to <- NULL

  # collect milestone network
  lineages <- slingshot::slingLineages(sds)
  dist <- igraph::distances(slingshot::slingMST(sds))
  cluster_network <- lineages |>
    purrr:::map_df(~ tibble::tibble(from = .[-length(.)], to = .[-1])) |>
    unique() |>
    dplyr:::mutate(
      length = dist[cbind(from, to)],
      directed = TRUE
    )

  # collect dimred
  dimred <- slingshot::slingReducedDim(sds)

  # collect progressions
  lin_assign <- apply(slingshot::slingCurveWeights(sds), 1, which.max)

  progressions <- purrr::map_df(seq_along(lineages), function(l) {
    ind <- lin_assign == l
    lin <- lineages[[l]]
    pst.full <- slingshot::slingPseudotime(sds, na = FALSE)[,l]
    pst <- pst.full[ind]
    means <- sapply(lin, function(clID){
      stats::weighted.mean(pst.full, cluster[,clID])
    })
    non_ends <- means[-c(1,length(means))]
    edgeID.l <- as.numeric(cut(pst, breaks = c(-Inf, non_ends, Inf)))
    from.l <- lineages[[l]][edgeID.l]
    to.l <- lineages[[l]][edgeID.l + 1]
    m.from <- means[from.l]
    m.to <- means[to.l]

    pct <- (pst - m.from) / (m.to - m.from)
    pct[pct < 0] <- 0
    pct[pct > 1] <- 1

    tibble::tibble(cell_id = names(which(ind)), from = from.l, to = to.l, percentage = pct)
  })

  output <-
    dynwrap::wrap_data(
      cell_ids = rownames(sds)) |>
    dynwrap::add_trajectory(
      milestone_network = cluster_network,
      progressions = progressions) |>
    dynwrap::add_dimred(
      dimred = dimred)

  output
}

clean_pickle = function(path) {
  olds = list.files(dirname(path), '_[0-9a-f]{32}[.]pickle$', full.names = TRUE)
  olds = c(olds, path)  # `path` may not exist; make sure it is in target paths
  base = basename(olds)
  keep = basename(path) == base  # keep this file (will cache to this file)
  base = substr(base, 1, nchar(base) - 37)  # 37 = 1 (_) + 32 (md5 sum) + 4 (.rds)
  unlink(olds[(base == base[keep][1]) & !keep])
}

cache_pickle <- function(
  expr = {}, rerun = FALSE,
  file = 'cache.rds', dir = 'cache/',
  hash = NULL, clean = TRUE,
  ...
) {
  if (xfun::loadable('knitr')) {
    if (missing(file) && !is.null(lab <- knitr::opts_current$get('label')))
      file = paste0(lab, '.rds')
    if (missing(dir) && !is.null(d <- knitr::opts_current$get('cache.path')))
      dir = d
  }

  path = paste0(dir, file)
  if (!grepl(r <- '([.]pickle)$', path)) path = paste0(path, '.pickle')

  code = deparse(substitute(expr))
  md5  = xfun:::md5sum_obj(code)

  if (identical(hash, 'auto')) hash = xfun:::global_vars(code, parent.frame(2))
  if (is.list(hash)) md5 = xfun:::md5sum_obj(c(md5, xfun:::md5sum_obj(hash)))

  path = sub(r, paste0('_', md5, '\\1'), path)
  if (rerun) unlink(path)
  if (clean) clean_pickle(path)

  library(reticulate)
  pickle <- import("pickle")
  builtins <- import_builtins()
  if (xfun::file_exists(path)) {
    with(builtins$open(path, "rb") %as% file,
         pickle$load(file))
  } else {
    obj = expr  # lazy evaluation
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    with(builtins$open(path, "wb") %as% file,
         pickle$dump(obj, file))
    obj
  }
}

readPickle <- function(path) {
  library(reticulate)
  pickle <- import("pickle")
  builtins <- import_builtins()
  with(builtins$open(path, "rb") %as% file,
       pickle$load(file))
}

#' Create Paga Object from AnnData
#'
#' @param adata \code{reticulate} reference to python AnnData object.
#' @param basis Dimensionality reduction name
createPagaObject <- function(adata, basis = "umap") {
  group_name <- adata$uns['paga']$groups
  obs <- adata$obs
  group_lvs <- levels(obs[[group_name]])

  paga <- list(
    connectivities = adata$uns['paga']$connectivities,
    connectivities_tree = adata$uns['paga']$connectivities_tree,
    transitions_confidence = adata$uns['paga']$transitions_confidence,
    threshold = adata$uns["paga"]$threshold,
    group_name = group_name,
    groups = group_lvs,
    group_colors = setNames(
      ggthemes::tableau_color_pal('Tableau 20')(length(group_lvs)), group_lvs),
    embeddings = data.frame(
      Dim_1 = as.data.frame(adata$obsm[paste0("X_", basis)])[, 1L],
      Dim_2 = as.data.frame(adata$obsm[paste0("X_", basis)])[, 2L],
      group = obs[[group_name]]
    )
  )

  for (mtx_name in c("connectivities", "connectivities_tree",
                     "transitions_confidence")) {
    if (!is.null(paga[[mtx_name]]))
      paga[[mtx_name]] <- as.matrix(paga[[mtx_name]])
  }

  # Calculate paga positions on cell embeddings
  paga$position <- dplyr::group_by(paga$embeddings, group) %>%
    dplyr::summarise(x = median(Dim_1), y = median(Dim_2))

  rownames(paga$connectivities) <- c(1:nrow(paga$position))
  colnames(paga$connectivities) <- c(1:nrow(paga$position))

  paga
}

#' Plot PAGA
#'
#' @param paga Paga object
#' @param threshold Do not draw edges for weights below this threshold.
#' Set to 0 if you want all edges. Discarding low-connectivity edges
#' helps in getting a much clearer picture of the graph.
plotPAGA <- function(paga, threshold = 0, edge_cex = 3, node_cex = 2, label_cex = 2) {
  # Inspired by https://romanhaa.github.io/blog/paga_to_r/
  paga_edges <- tibble::tibble(
    group1 = rownames(paga$connectivities)[row(paga$connectivities)[upper.tri(paga$connectivities)]],
    group2 = colnames(paga$connectivities)[col(paga$connectivities)[upper.tri(paga$connectivities)]],
    weight = paga$connectivities[upper.tri(paga$connectivities)]
  ) %>%
    dplyr::mutate(
      x1 = paga$position$x[match(.$group1, rownames(paga$position))],
      y1 = paga$position$y[match(.$group1, rownames(paga$position))],
      x2 = paga$position$x[match(.$group2, rownames(paga$position))],
      y2 = paga$position$y[match(.$group2, rownames(paga$position))]
    )

  paga_edges_sig <- dplyr::filter(paga_edges, weight >= threshold)
  paga_edges_und <- dplyr::filter(paga_edges, weight < threshold)

  p <- ggplot2::ggplot(
      paga$position, ggplot2::aes(x, y)) +
    ggplot2::geom_point(
      data = paga$embeddings,
      mapping = ggplot2::aes(
        x = Dim_1, y = Dim_2, color = group),
      shape = 19,
      alpha = 0.1,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = paga_edges_und,
      mapping = ggplot2::aes(
        x = x1, y = y1, xend = x2, yend = y2),
      linetype = "dashed",
      color = "grey",
      size = paga_edges_und$weight*edge_cex,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = paga_edges_sig,
      mapping = ggplot2::aes(
        x = x1, y = y1, xend = x2, yend = y2),
      linetype = "solid",
      color = "black",
      size = paga_edges_sig$weight*edge_cex,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = group),
      size = 3*node_cex, alpha = 1, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = paga$group_colors) +
    ggplot2::geom_text(
      ggplot2::aes(label = group), color = "black",
      size = 3*label_cex, fontface = "bold") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )

  p
}

plotPagaGraph <- function(adata, basis = "umap", threshold = 0, ...) {
  paga <- createPagaObject(adata, basis = basis)
  plotPAGA(paga, threshold = threshold, ...)
}

plotPagaCompare <- function(adata, basis = "umap", threshold = 0) {
  paga <- createPagaObject(adata, basis = basis)

  obs <- adata$obs

  embed_df <- tibble::tibble(
    Dim_1 = paga$embeddings$Dim_1,
    Dim_2 = paga$embeddings$Dim_2,
    group = obs[[paga$group_name]]
  )

  group_centers <- embed_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(x = median(Dim_1), y = median(Dim_2))

  embed_plot <- embed_df %>%
    ggplot2::ggplot(ggplot2::aes(Dim_1, Dim_2, color = group)) +
    ggplot2::geom_point(size = 0.1, show.legend = FALSE) +
    ggplot2::geom_text(
      data = group_centers,
      mapping = ggplot2::aes(x, y, label = group),
      size = 4.5,
      color = "black",
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = paga$group_colors) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )

  paga_plot <- plotPAGA(paga, threshold = threshold)

  patchwork::wrap_plots(embed_plot, paga_plot)
}

plotPagaArrow <- function(
    adata, basis = "umap", threshold = 0,
    edge_cex = 3, node_cex = 2, label_cex = 2,
    arrow_gap = 0.5, point_alpha = 0.1) {
  paga <- createPagaObject(adata, basis = basis)

  paga_edges <- tibble::tibble(
    group1 = rownames(paga$connectivities)[row(paga$connectivities)[upper.tri(paga$connectivities)]],
    group2 = colnames(paga$connectivities)[col(paga$connectivities)[upper.tri(paga$connectivities)]],
    weight = paga$connectivities[upper.tri(paga$connectivities)]
  ) %>%
    dplyr::mutate(
      x1 = paga$position$x[match(.$group1, rownames(paga$position))],
      y1 = paga$position$y[match(.$group1, rownames(paga$position))],
      x2 = paga$position$x[match(.$group2, rownames(paga$position))],
      y2 = paga$position$y[match(.$group2, rownames(paga$position))]
    )

  paga_arrows <- tibble::tibble(
    group1 = character(),
    group2 = character(),
    weight = double(),
    x1 = double(),
    y1 = double(),
    x2 = double(),
    y2 = double()
  )
  
  paga_arrows_dimensions <- nrow(paga$transitions_confidence)

  comparisons_done <- c()

  for ( i in 1:paga_arrows_dimensions ) {
    # loop through columns
    for ( j in 1:paga_arrows_dimensions ) {
      # skip cell if on diagonal
      if ( i == j ) {
        next
        # skip cell if transition between these clusters has already been extracted
      } else if ( paste0(i, "/", j) %in% comparisons_done | paste0(j, "/", i) %in% comparisons_done ) {
        next
        # if none of the above, go ahead
      } else {
        # get value for transition i to j
        i_to_j <- paga$transitions_confidence[j,i]
        # get value for transition j to i (other side of diagonal)
        j_to_i <- paga$transitions_confidence[i,j]
        # if i to j is more confident than j to i
        if ( i_to_j > j_to_i ) {
          x1 <- paga$position$x[i]
          y1 <- paga$position$y[i]
          x2 <- paga$position$x[j]
          y2 <- paga$position$y[j]
          x_length = x2 - x1
          y_length = y2 - y1
          gap = arrow_gap / sqrt(x_length ^ 2 + y_length ^ 2)
          x1_new = x1 + gap * x_length
          y1_new = y1 + gap * y_length
          x2_new = x1 + (1 - gap) * x_length
          y2_new = y1 + (1 - gap) * y_length
          new_entry <- tibble::tibble(
            group1 = (i - 1),
            group2 = (j - 1),
            weight = i_to_j,
            x1 = x1_new,
            y1 = y1_new,
            x2 = x2_new,
            y2 = y2_new
          )
          paga_arrows <- rbind(
            paga_arrows,
            new_entry
          )
          # if j to i is more confident than i to j
        } else if ( j_to_i > i_to_j ) {
          x1 <- paga$position$x[j]
          y1 <- paga$position$y[j]
          x2 <- paga$position$x[i]
          y2 <- paga$position$y[i]
          x_length = x2 - x1
          y_length = y2 - y1
          gap = arrow_gap / sqrt(x_length ^ 2 + y_length ^ 2)
          x1_new = x1 + gap * x_length
          y1_new = y1 + gap * y_length
          x2_new = x1 + (1 - gap) * x_length
          y2_new = y1 + (1 - gap) * y_length
          new_entry <- tibble::tibble(
            group1 = (j - 1),
            group2 = (i - 1),
            weight = j_to_i,
            x1 = x1_new,
            y1 = y1_new,
            x2 = x2_new,
            y2 = y2_new
          )
          paga_arrows <- rbind(
            paga_arrows,
            new_entry
          )
          # else (for example both 0)
        } else {
          next
        }
        # add comparison to list of already performed comparisons
        comparisons_done <- c(comparisons_done, paste0(i,"/",j))
      }
    }
  }

  paga_edges <- paga_edges %>%
    dplyr::filter(weight >= threshold)
  
  paga_arrows <- paga_arrows %>%
    dplyr::filter(weight >= threshold)
  
  plot <- ggplot2::ggplot()

  plot <- plot + ggplot2::geom_point(
    data = paga$embeddings,
    mapping = ggplot2::aes(
      x = Dim_1, y = Dim_2, color = group),
    shape = 19,
    alpha = point_alpha,
    show.legend = FALSE
  )

  # Plot grey edges
  plot <- plot + ggplot2::geom_segment(
    data = paga_edges,
    mapping = ggplot2::aes(
      x = x1, y = y1, xend = x2, yend = y2),
    size = paga_edges$weight*edge_cex,
    linetype = "dashed",
    colour = "grey",
    alpha = 0.75,
    show.legend = FALSE
  )

  # Plot paga node
  plot <- plot + ggplot2::geom_point(
    data = paga$position,
    mapping = ggplot2::aes(x, y, color = group),
    size = 3*edge_cex,
    show.legend = FALSE
  )

  # Plot arrow
  plot <- plot + ggplot2::geom_segment(
    data = paga_arrows,
    mapping = ggplot2::aes(
      x = x1, y = y1, xend = x2, yend = y2),
    size = paga_arrows$weight*edge_cex,
    arrow = ggplot2::arrow(
      length = ggplot2::unit(0.004*edge_cex, "npc"),
      type = "closed", angle = 15),
    show.legend = FALSE
  )

  # Plot node label
  plot <- plot + ggplot2::geom_text(
    data = paga$position,
    mapping = ggplot2::aes(x, y, label = group),
    color = "black",
    fontface = "bold",
    size = 3*label_cex
  )

  plot <- plot + ggplot2::scale_color_manual(
      values = paga$group_colors) +
    ggplot2::labs(x = "Dim_1", y = "Dim_2") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )

  plot
}
