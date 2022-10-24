suppressMessages(library(SingleCellExperiment))

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

runUniTVelo <- function(sce, dimred = NULL, use.dimred = NULL, ...) {
  if (is.null(dimred)) {
    if (!is.null(use.dimred)) {
      dimred <- reducedDim(x, use.dimred)
    }
  }

  spliced <- assay(sce, "spliced")
  unspliced <- assay(sce, "unspliced")
  X <- assay(sce, "counts")
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

wrapSlingshotToDynverse <- function(sds) {
  cluster <- slingshot::slingClusterLabels(sds)
  start_cell <- apply(slingshot::slingPseudotime(sds), 1, min) %>%
    sort() %>% head(1) %>% names()
  start.clus <- names(which(cluster[start_cell,] == 1))

  # satisfy r cmd check
  from <- to <- NULL

  # collect milestone network
  lineages <- slingshot::slingLineages(sds)
  dist <- igraph::distances(slingshot::slingMST(sds))
  cluster_network <- lineages %>%
    purrr:::map_df(~ tibble::tibble(from = .[-length(.)], to = .[-1])) %>%
    unique() %>%
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
      cell_ids = rownames(sds)) %>%
    dynwrap::add_trajectory(
      milestone_network = cluster_network,
      progressions = progressions) %>%
    dynwrap::add_dimred(
      dimred = dimred)

  output
}
