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

plot_markers <- function(obj, ...) {
  Seurat::DotPlot(
    obj,
    features = list(
      # Epidermal Cells
      "FDH" = c(
        "TraesCS4B02G297500",
        "TraesCS4D02G296400",
        "TraesCS4A02G007400"
      ),
      "ATML1" = c(
        "TraesCS2A02G474000",
        "TraesCS2D02G473700"
      ),
      "DCR" = c(
        "TraesCS1A02G341300",
        "TraesCS1D02G343400"
      ),
      
      # EP3 是排水孔相关基因
      "EP3" = c(
        "TraesCS2A02G350700",
        "TraesCS2D02G348800",
        "TraesCS6D02G199500",
        "TraesCS6A02G216100"
      ),
      
      # Guardian Cells
      "ALMT12" = c(
        "TraesCS1D02G194000",
        "TraesCS1A02G189900",
        "TraesCS1B02G192000"
      ),
      "MYB60" = c(
        "TraesCS4A02G322200",
        "TraesCS5D02G552200"
      ),
      "HIC" = c(
        "TraesCS4D02G226100"
      ),
      
      # Mesophyll Cells
      "RBCS" = c(
        "TraesCS2A02G066800",
        #"TraesCS2B02G079200",
        #"TraesCS2D02G065200",
        #"TraesCS2D02G065300",
        "TraesCS5A02G165400",
        #"TraesCS5A02G165700",
        #"TraesCS5B02G162600",
        #"TraesCS5B02G162800",
        #"TraesCS5D02G169600",
        "TraesCS5D02G169900"
      ),
      "CAB3" = c(
        "TraesCS7A02G276400",
        #"TraesCS1D02G411300",
        "TraesCS1B02G317500",
        #"TraesCS7D02G276300",
        #"TraesCS5B02G353200",
        #"TraesCS5A02G350600",
        "TraesCS1A02G403300"
      ),
      "LHCB2.1" = c(
        "TraesCS5D02G329200",
        "TraesCS5B02G322900",
        "TraesCS5A02G322500"
      ),
      "CA1" = c(
        #"TraesCS7D02G443400",
        "TraesCS7B02G354800",
        "TraesCS3A02G230000",
        #"TraesCS3D02G223300",
        "TraesCS3B02G259300"
      ),
      "AOC2" = c(
        "TraesCS6D02G314300",
        "TraesCS6A02G334800",
        "TraesCS6B02G365200"
      ),
      
      # Vascular Cells
      "SULTR3;4" = c(
        "TraesCS7A02G088700",
        "TraesCS4A02G388000",
        "TraesCS7D02G084100"
      ),
      "TaGSr" = c(
        "TraesCS4B02G240900",
        "TraesCS4D02G240700",
        "TraesCS4A02G063800"
      ),
      "gl-OXO" = c(
        "TraesCS4D02G032000",
        "TraesCS4B02G033300",
        "TraesCS4A02G279200",
        "TraesCS4D02G031800"
      ),
      "TaSUT1" = c(
        "TraesCS4A02G016400",
        "TraesCS4B02G287800",
        "TraesCS4D02G286500"
      ),
      "CPIII" = c(
        "TraesCS6B02G050700",
        "TraesCS6D02G041700",
        "TraesCS6A02G036100"
      ),
      
      # Cortex Cells
      "AT1G62510" = c(
        "TraesCS2A02G424800",
        "TraesCS2B02G444500",
        #"TraesCS2D02G422700",
        #"TraesCS2D02G422800",
        "TraesCS2A02G424861"
      ),
      # Bulliform
      "PFA-DSP2" = c(
        "TraesCS5B02G163200",
        "TraesCS5D02G170400"
      )
    ),
    ...
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text.x = ggplot2::element_text(angle = 90)
    )
}

choose_indent_res <- function(obj, meta_col) {
  Seurat::Idents(obj) <- meta_col
  lvs <- type.convert(levels(obj))
  Seurat::Idents(obj) <- factor(Seurat::Idents(obj), levels = sort(lvs))
  obj
}