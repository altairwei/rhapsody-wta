library(magrittr)

#' Perform ORA enrichment analysis for DS results
#'
#' @param ds DS results
#' @param contrasts Which contrast to pull
#' @param clusters Which clusters to compare
#' @param ont Gene ontology to display
#' @param simplify Simplify GO term
#'
#' @return A list of two elements: \code{enr} and \code{fc}
#'
performDSEnrichORA <- function(
  ds, contrasts, clusters,
  ont = c("BP", "MF", "CC", "ALL"),
  simplify = TRUE
) {
  ont <- match.arg(ont)
  
  # We use all detected genes in DS analysis of given cluster as background
  background <- rhapsodykit::diff_state_pull(
    ds, contrasts, clusters, "gene")
  
  sig_genes <- ds %>%
    rhapsodykit::diff_state_significant()
  genes <- sig_genes %>%
    rhapsodykit::diff_state_pull(contrasts, clusters, "gene")
  fc <- sig_genes %>%
    rhapsodykit::diff_state_pull(contrasts, clusters, "logFC")
  names(fc) <- genes
  
  enr <- clusterProfiler::enrichGO(
    genes,
    org.Taestivum.iwgsc.db,
    keyType = "GID",
    ont = ont,
    universe = background
  )
  
  if (simplify && ont != "ALL")
    enr <- clusterProfiler::simplify(enr, cutoff = 0.5)
  
  list(
    enr = enr,
    fc = fc
  )
}

#' Compare ORA enrichment results
#'
#' @inheritParams performDSEnrichORA
#'
#' @return A list of two elements: \code{enr} and \code{fc}
#'
performDSCompareORA <- function(
  ds, contrasts, clusters,
  ont = c("BP", "MF", "CC", "ALL"),
  simplify = TRUE
) {
  ont <- match.arg(ont)
  
  if (is.null(names(contrasts)))
    names(contrasts) <- contrasts
  if (is.null(names(clusters)))
    names(clusters) <- clusters
  
  background <- purrr::imap(clusters, function(clr, clr_name) {
    purrr::imap(contrasts, function(con, con_name) {
      ds %>% rhapsodykit::diff_state_pull(con, clr, "gene")
    })
  }) %>%
    unlist() %>%
    unique()
  
  compare_df <- purrr::imap(clusters, function(clr, clr_name) {
    purrr::imap(contrasts, function(con, con_name) {
      df <- ds %>%
        rhapsodykit::diff_state_significant() %>%
        rhapsodykit::diff_state_pull(con, clr, c("gene", "logFC"))
      df$gene
      data.frame(
        gene = df$gene,
        contrast = con_name,
        cluster = clr_name
      )
    })
  }) %>%
    unlist(recursive = FALSE) %>%
    dplyr::bind_rows()
  
  compare_df$cluster <- factor(compare_df$cluster, levels = names(clusters))
  compare_df$contrast <- factor(compare_df$contrast, levels = names(contrasts))
  
  comp_enr <- clusterProfiler::compareCluster(
    gene ~ cluster + contrast,
    data = compare_df,
    fun = "enrichGO",
    OrgDb = org.Taestivum.iwgsc.db,
    keyType = "GID",
    ont = ont,
    universe = background
  )
  
  if (simplify && ont != "ALL")
    comp_enr <- clusterProfiler::simplify(comp_enr, cutoff = 0.5)
  
  comp_enr
}

performDSEnrichGSEA <- function(
  ds, contrasts, clusters,
  ont = c("BP", "MF", "CC", "ALL"),
  simplify = TRUE
) {
  ont <- match.arg(ont)
  
  ## 我们不需要依据 FDR 来筛选基因，只需要 logFC 来排序
  gsea_df <- ds %>%
    rhapsodykit::diff_state_pull(contrasts, clusters, c("gene", "logFC"))
  
  ## feature 1: numeric vector
  geneList <- gsea_df[, 2]
  
  ## feature 2: named vector
  names(geneList) <- as.character(gsea_df[, 1])
  
  ## feature 3: decreasing order
  geneList <- sort(geneList, decreasing = TRUE)
  
  gsea <- clusterProfiler::gseGO(
    geneList,
    ont = ont,
    OrgDb = org.Taestivum.iwgsc.db,
    keyType = "GID"
  )
  
  if (simplify && ont != "ALL")
    gsea <- clusterProfiler::simplify(gsea)
  
  list(
    gsea = gsea,
    fc = geneList
  )
  
}

performDSCompareGSEA <- function(
  ds, contrasts, clusters,
  ont = c("BP", "MF", "CC", "ALL"),
  simplify = TRUE
) {
  ont <- match.arg(ont)
  
  if (is.null(names(contrasts)))
    names(contrasts) <- contrasts
  if (is.null(names(clusters)))
    names(clusters) <- clusters
  
  compare_df <- purrr::imap(clusters, function(clr, clr_name) {
    purrr::imap(contrasts, function(con, con_name) {
      df <- ds %>%
        rhapsodykit::diff_state_pull(con, clr, c("gene", "logFC"))
      geneList <- df[, 2]
      names(geneList) <- as.character(df[, 1])
      geneList <- sort(geneList, decreasing = TRUE)
      data.frame(
        gene = names(geneList),
        fc = geneList,
        cluster = clr_name,
        contrast = con_name
      )
    })
  }) %>%
    unlist(recursive = FALSE) %>%
    dplyr::bind_rows()
  
  compare_df$cluster <- factor(compare_df$cluster, levels = names(clusters))
  compare_df$contrast <- factor(compare_df$contrast, levels = names(contrasts))
  
  comp_enr <- clusterProfiler::compareCluster(
    gene | fc ~ cluster + contrast,
    data = compare_df,
    OrgDb = org.Taestivum.iwgsc.db,
    fun = "gseGO",
    keyType = "GID",
    ont = ont
  )
  
  if (simplify && ont != "ALL")
    comp_enr <- clusterProfiler::simplify(comp_enr, cutoff = 0.5)
  
  comp_enr
}

#' Perform comparative ORA for cell type marker genes
#'
#' @inheritParams performDSEnrichORA
#' @param markerlist A list of character vectors.
#' @param clusters A (named) character vector which represents
#' selected cell types.
#' @param background Background genes used for ORA.
#' @return Results of \code{clusterProfiler::compareCluster}
#'
performMarkerCompareORA <- function(
    markerlist, clusters, background,
    ont = c("BP", "MF", "CC", "ALL"),
    simplify = TRUE
) {
  ont <- match.arg(ont)
  
  if (is.null(names(clusters)))
    names(clusters) <- clusters
  
  compare_df <- purrr::imap(clusters, function(clr, clr_name) {
    df <- data.frame(
      gene = markerlist[[clr]],
      cluster = rep(clr_name, length(markerlist[[clr]]))
    )
    df
  }) %>%
    dplyr::bind_rows()
  
  compare_df$cluster <- factor(compare_df$cluster, levels = names(clusters))
  
  comp_enr <- clusterProfiler::compareCluster(
    gene ~ cluster,
    data = compare_df,
    fun = "enrichGO",
    OrgDb = org.Taestivum.iwgsc.db,
    keyType = "GID",
    ont = ont,
    universe = background
  )
  
  if (simplify && ont != "ALL")
    comp_enr <- clusterProfiler::simplify(comp_enr, cutoff = 0.5)
  
  comp_enr
}

#' Print comparative ORA enrichment results as table
#'
#' @param enr Comparative ORA enrichment results.
#' @param heading Markdown heading symbols.
print_enrich_table <- function(enr, heading) {
  df <- enr@compareClusterResult |>
    dplyr::select(
      cluster, ID, Description,
      Count, p.adjust, geneID) |>
    dplyr::group_by(cluster)
  
  df_list <- dplyr::group_split(df)
  names(df_list) <- dplyr::group_keys(df)$cluster
  
  for (clr in names(df_list)) {
    cat(heading, clr, "\n\n")
    print(
      htmltools::tagList(
        reactable::reactable(
          df_list[[clr]],
          columns = list(
            "cluster" = reactable::colDef(maxWidth = 80),
            "ID" = reactable::colDef(maxWidth = 120),
            "Description" = reactable::colDef(minWidth = 100),
            "Count" = reactable::colDef(maxWidth = 80),
            "p.adjust" = reactable::colDef(
              maxWidth = 80, format = reactable::colFormat(digits = 2)),
            "geneID" = reactable::colDef(minWidth = 250)
          )
        )
      )
    )
    cat("\n\n")
  }
}

performDACompareORA <- function(da, sets,
                            ont = c("BP", "MF", "CC", "ALL"),
                            simplify = TRUE,
                            logFC = c("both", "up", "down")
) {
  stopifnot(is.list(sets) && !is.null(names(sets)) && length(sets) > 1)
  logFC = match.arg(logFC)

  compare_list <- sets %>%
    lapply(function(set) da[[ set[[1]] ]][[ set[[2]] ]])
  
  comp_enr <- clusterProfiler::compareCluster(
    compare_list, fun = "enrichGO",
    OrgDb = org.Taestivum.iwgsc.db,
    keyType = "GID",
    ont = ont)

  if (simplify && ont != "ALL")
    comp_enr <- clusterProfiler::simplify(comp_enr)

  comp_enr
}