#' Extract information from SingleCellExperiment object to create CNA object
#'
#' @param x SingleCellExperiment object.
#' @param samplem_key String denoting the name of the sample-level identifier (e.g. DonorID). 
#' @param samplem_vars Which sample-level covariates to include. 
#' @param dimred Dimensionality reduction to use as input for building the NN graph.
#'
create_object.SingleCellExperiment <- function(x, samplem_key, samplem_vars, dimred = "PCA", ...) {
  samplem_vars <- c(samplem_vars, samplem_key)
  colDf <- as.data.frame(colData(x))

  samplem_df <- colDf |>
    dplyr::select(tidyselect::one_of(samplem_vars)) |>
    unique() |>
    tibble::remove_rownames()

  obs_df <- tibble::rownames_to_column(colDf, 'CellID')
  
  if (nrow(samplem_df) == nrow(obs_df)) {
    stop(
      'Sample-level metadata is same length as cell-level metadata.       
         Please check that samplem_vars are sample-level covariates.'
    )
  }

  nngraph <- Seurat::FindNeighbors(
    reducedDim(x, dimred)[, 1:10], compute.SNN = FALSE,
    return.neighbor = FALSE, ...)$nn

  rcna_object <- list(
    samplem = samplem_df,
    obs = obs_df, 
    connectivities = nngraph,
    samplem_key = samplem_key,
    obs_key = 'CellID',
    N = nrow(samplem_df)
  )

  rcna_object
}

#' Perform co-varing neighborhood analysis on SingleCellExperiment object
#'
#' @param x SingleCellExperiment object.
#' @param test_var Contrast variable to test for association. 
#' @param samplem_key String denoting the name of the sample-level identifier (e.g. DonorID). 
#' @param dimred Dimensionality reduction to use as input for building the NN graph.
#' @param batches Name of batch variable. Currently only one categorical variable allowed. 
#' @param covs Name(s) of other (numerical) covariates to control for. 
#' @param nsteps TBD
#' @param verbose TBD
#'
association.SingleCellExperiment <- function(
    x, test_var, samplem_key, dimred = "PCA", 
    batches = NULL, covs = NULL, nsteps = NULL, verbose=TRUE, 
    assay = NULL) {
  
  ## (1) format data 
  covs_keep <- test_var
  if (!is.null(batches)) covs_keep <- c(covs_keep, batches)
  if (!is.null(covs)) covs_keep <- c(covs_keep, covs)
  
  rcna_data <- create_object.SingleCellExperiment(
    x = x,
    samplem_key = samplem_key,
    samplem_vars = covs_keep,
    dimred = dimred)

  yvals <- rcna_data$samplem[[test_var]]
  if (is(yvals, 'character') | is(yvals, 'factor') | is(yvals, 'integer') ) {
    stop('test_var is of class {class(yvals)}. It must be numeric variable for association testing.')
  }

  ## (2) do association
  cna_res <- rcna::association(
    data = rcna_data, 
    y = yvals,
    batches = batches,
    covs = covs,
    nsteps = nsteps,
    suffix = '',
    force_recompute = TRUE,
    return_nam = TRUE,
    verbose = verbose
  )
  
  ## (3) save results
  metadata(x)$cnaRes <- cna_res
  reducedDim(x, "CNA") <- cna_res$NAM_embeddings
  colnames(reducedDim(x, "CNA")) <- NULL

  x$cna_ncorrs <- cna_res$ncorrs[colnames(x), , drop=TRUE]
  x$cna_ncorrs_fdr05 <- rep(0, ncol(x))
  if (!is.null(cna_res$fdr_5p_t)) {
    idx_passed <- which(abs(x$cna_ncorrs) >= cna_res$fdr_5p_t)
    x$cna_ncorrs_fdr05[idx_passed] <- x$cna_ncorrs[idx_passed]
  }

  x$cna_ncorrs_fdr10 <- rep(0, ncol(x))
  if (!is.null(cna_res$fdr_10p_t)) {
    idx_passed <- which(abs(x$cna_ncorrs) >= cna_res$fdr_10p_t)
    x$cna_ncorrs_fdr10[idx_passed] <- x$cna_ncorrs[idx_passed]
  }

  x
}

plotNcorrRidges <- function(
    sce, theme_size = 11,
    y_expand = ggplot2::expansion(mult = c(0, 0.2)),
    group_by = "cluster_id") {
  cnaRes <- metadata(sce)$cnaRes
  data.frame(cna_ncorrs = sce$cna_ncorrs,
             cellType = colData(sce)[, group_by]) |>
    ggplot2::ggplot(ggplot2::aes(
      x = cna_ncorrs, y = forcats::fct_reorder(cellType, cna_ncorrs),
      fill = cellType)) +
    ggplot2::ylab("cellType") +
    ggridges::geom_density_ridges(scale = 4) +
    ggplot2::geom_vline(xintercept = cnaRes$fdr_5p_t, linetype="dashed") +
    ggplot2::geom_vline(xintercept = -cnaRes$fdr_5p_t, linetype="dashed") +
    ggplot2::annotate(
      geom = "label", x = cnaRes$fdr_5p_t, y = Inf, 
      label = "5% FDR", hjust = 0.5, vjust = 1, size = 5) +
    ggplot2::annotate(
      geom = "label", x = -cnaRes$fdr_5p_t, y = Inf, 
      label = "5% FDR", hjust = 0.5, vjust = 1, size = 5) +
    ggplot2::scale_y_discrete(limits=rev, expand = y_expand) +
    ggridges::theme_ridges(theme_size)
}

coffGeneList <- function(ncorrs, exprs) {
  ncorrs <- as(matrix(ncorrs, ncol = 1), "dgCMatrix")
  genes <- qlcMatrix::corSparse(ncorrs, exprs)
  genes <- structure(as.numeric(genes), names = colnames(exprs))
  genes <- sort(genes, decreasing = TRUE)
  genes
}

runGSEA <- function(x) {
  rankLists <- x
  stopifnot(length(rankLists) > 0)
  
  if (length(rankLists) == 1) {
    gsea <- clusterProfiler::gseGO(
      geneList = rankLists[[1]],
      ont = "BP",
      OrgDb = 'org.Taestivum.iwgsc.db',
      keyType = "GID"
    )
  } else {
    gsea <- clusterProfiler::compareCluster(
      geneClusters = rankLists,
      fun = "gseGO",
      ont = "BP",
      OrgDb = 'org.Taestivum.iwgsc.db',
      keyType = "GID"
    )
  }
  
  gsea <- clusterProfiler::simplify(gsea)
  
  gsea
}


gseaTable <- function(gsea) {
  digits_2 <- reactable::colFormat(digits = 2)
  gsea@compareClusterResult %>%
    dplyr::select(Cluster, ID, Description, setSize,
                  enrichmentScore, NES, pvalue, p.adjust,
                  qvalues, rank, leading_edge) %>%
    reactable::reactable(columns = list(
      "ID" = reactable::colDef(minWidth = 120),
      "Description" = reactable::colDef(minWidth = 200),
      "enrichmentScore" = reactable::colDef(format = digits_2, maxWidth = 80),
      "NES" = reactable::colDef(format = digits_2, maxWidth = 70),
      "pvalue" = reactable::colDef(format = digits_2, maxWidth = 50),
      "p.adjust" = reactable::colDef(format = digits_2, maxWidth = 50),
      "qvalues" = reactable::colDef(format = digits_2, maxWidth = 50)
    ))
}

plot_bootstrap_distribution_for_mock <- function(
    res, clusters = NULL,
    facet_by = "cellTypes", ncol = NULL, nrow = NULL,
    theme_size = 11
) {
  # Use BCa as CI
  df_to_plot <- res$results %>%
    tibble::as_tibble() %>%
    dplyr::filter(method == "BCa") %>%
    dplyr::mutate(
      time = sapply(strsplit(cond, split = "-"), "[", 1),
      rep = sapply(strsplit(as.character(subject), split = "-"),
                   function(x) paste(x[2:3], collapse = "-")),
      treatment = sapply(strsplit(cond, split = "-"), "[", 2)
    )
  
  prop_df <- res$thetastar %>%
    as.data.frame()
  
  colnames(prop_df) <- paste("BS_", seq_len(ncol(res$thetastar)), sep = "")
  
  df_to_plot <- dplyr::bind_cols(df_to_plot, prop_df) %>%
    tidyr::pivot_longer(tidyr::starts_with("BS_"), names_to = "bootstrap", values_to = "prop")
  
  if (!is.null(clusters)) {
    df_to_plot <- dplyr::filter(df_to_plot, cellTypes %in% clusters)
    df_to_plot$cellTypes <- factor(df_to_plot$cellTypes, levels = clusters)
  }
  
  p <- ggplot2::ggplot(df_to_plot,
                       ggplot2::aes(x = time, y = prop, fill = cond)) +
    ggplot2::geom_boxplot(
      width = 0.2, outlier.shape = NA,
      position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(group = treatment),
      stat = "summary", fun = median,
      position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_bw(theme_size) +
    NULL
  
  if (!is.null(facet_by))
    p <- p + ggplot2::facet_wrap(facet_by, scales = "free_y", ncol = ncol, nrow = nrow, drop = FALSE)
  
  p
}

plot_bootstrap_distribution <- function(
    res, clusters = NULL,
    facet_by = "cellTypes", ncol = NULL, nrow = NULL
) {
  # Use BCa as CI
  df_to_plot <- res$results %>%
    tibble::as_tibble() %>%
    dplyr::filter(method == "BCa") %>%
    dplyr::mutate(
      time = sapply(strsplit(cond, split = "-"), "[", 1),
      rep = sapply(strsplit(as.character(subject), split = "-"),
                   function(x) paste(x[2:3], collapse = "-")),
      treatment = sapply(strsplit(cond, split = "-"), "[", 2)
    )
  
  prop_df <- res$thetastar %>%
    as.data.frame()
  
  colnames(prop_df) <- paste("BS_", seq_len(ncol(res$thetastar)), sep = "")
  
  df_to_plot <- dplyr::bind_cols(df_to_plot, prop_df) %>%
    tidyr::pivot_longer(tidyr::starts_with("BS_"), names_to = "bootstrap", values_to = "prop")
  
  #df_to_plot <- df_to_plot %>%
  #  dplyr::group_by(time, treatment, cellTypes) %>%
  #  dplyr::summarise(mean = mean(prop), sd = sd(prop))
  
  df_to_plot <- dplyr::bind_rows(
    df_to_plot,
    df_to_plot %>%
      dplyr::filter(time == "0DPI") %>%
      dplyr::mutate(treatment = "PNR2"),
    df_to_plot %>%
      dplyr::filter(time == "0DPI") %>%
      dplyr::mutate(treatment = "TR4")
  )
  
  if (!is.null(clusters)) {
    df_to_plot <- dplyr::filter(df_to_plot, cellTypes %in% clusters)
    df_to_plot$cellTypes <- factor(df_to_plot$cellTypes, levels = clusters)
  }
  
  p <- ggplot2::ggplot(df_to_plot,
                       ggplot2::aes(x = time, y = prop, color = treatment)) +
    ggplot2::geom_boxplot(
      width = 0.2, outlier.shape = NA,
      position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(group = treatment),
      stat = "summary", fun = median,
      position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_bw() +
    NULL
  
  if (!is.null(facet_by))
    p <- p + ggplot2::facet_wrap(facet_by, scales = "free_y", ncol = ncol, nrow = nrow, drop = FALSE)
  
  p
}
