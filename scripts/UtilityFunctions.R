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

#' Display cluster markers in a html table
#'
#' @param df Data frame returned by Seurat::FindAllMarkers
marker_table <- function(df) {
  df %>%
    dplyr::select(gene, dplyr::everything()) %>%
    reactable::reactable(
      defaultColDef = reactable::colDef(
        minWidth = 40
      ),
      columns = list(
        gene = reactable::colDef(minWidth = 80),
        p_val = reactable::colDef(
          cell = function(x) format(x, digits=3, scientific = TRUE)
        ),
        p_val_adj = reactable::colDef(
          cell = function(x) format(x, digits=3, scientific = TRUE)
        ),
        avg_log2FC = reactable::colDef(
          format = reactable::colFormat(
            digits = 2
          )
        ),
        pct.1 = reactable::colDef(
          format = reactable::colFormat(
            percent = TRUE,
            digits = 1
          )
        ),
        pct.2 = reactable::colDef(
          format = reactable::colFormat(
            percent = TRUE,
            digits = 1
          )
        ),
        .rownames = reactable::colDef(
          show = FALSE
        )
      )
    )
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