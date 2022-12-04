createCNAObject <- function(x, samplem_key, samplem_vars, dimred = "PCA", ...) {
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

