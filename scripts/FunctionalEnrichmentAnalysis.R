library(magrittr)

perform_ora <- function(ds, go_data, contrast, cluster) {
  df <- ds %>%
    rhapsodykit::diff_state_significant() %>%
    rhapsodykit::diff_state_pull(contrast, cluster, c("gene", "logFC"))

  go_data <- go_data %>%
    dplyr::filter(GO_ID != "")

  go2gene <- go_data %>%
    dplyr::select(GO_ID, Gene_ID)
  go2name <- go_data %>%
    dplyr::select(GO_ID, GO_Name)

  ora <- clusterProfiler::enricher(
    df$gene, TERM2GENE = go2gene, TERM2NAME = go2name)

  ora
}

prepare_significant_logFC <- function(ds, contrast, cluster) {
  df <- ds %>%
      rhapsodykit::diff_state_significant() %>%
      rhapsodykit::diff_state_pull(contrast, cluster, c("gene", "logFC"))

  fc_list <- df[, "logFC"]
  names(fc_list) <- as.character(df[, "gene"])

  fc_list
}

prepare_ranked_logFC <- function(ds, contrast, cluster) {
  df <-  ds %>%
    # We don't need to filter genes with given cutoff
    rhapsodykit::diff_state_pull(contrast, cluster, c("gene", "logFC"))

  fc_list <- df[, "logFC"]
  names(fc_list) <- as.character(df[, "gene"])

  ranked_list <- sort(fc_list, decreasing = TRUE)

  ranked_list
}

perform_gsea <- function(ds, go_data, contrast, cluster) {
  ranked_list <- prepare_ranked_logFC(ds, contrast, cluster)

  go_data <- go_data %>%
    dplyr::filter(GO_ID != "")

  go2gene <- go_data %>%
    dplyr::select(GO_ID, Gene_ID)
  go2name <- go_data %>%
    dplyr::select(GO_ID, GO_Name)

  gsea <- clusterProfiler::GSEA(
    ranked_list, TERM2GENE = go2gene, TERM2NAME = go2name)

  gsea
}

# Over Representation Analysis
# --------------------------------------

ora <- perform_ora(
  ds_res$DESeq2, go_data,
  "X1DPI.PNR2-X1DPI.MOCK",
  "Vascular Cells 1"
)

fc <- prepare_significant_logFC(
  ds_res$DESeq2, "X1DPI.PNR2-X1DPI.MOCK", "Vascular Cells 1")

ora %>%
  rhapsodykit::enrich_cnetplot(fold_change = fc, show_category = 20)

ora %>%
  enrichplot::pairwise_termsim() %>%
  enrichplot::emapplot(cex_label_category = 0.6)

# Gene Set Enrichment Analysis
# --------------------------------------

gsea <- perform_gsea(
  ds_res$DESeq2, go_data,
  "X1DPI.PNR2-X1DPI.MOCK",
  "Vascular Cells 1"
)

gsea %>%
  rhapsodykit::enrich_cnetplot(
    fold_change = gsea@geneList, show_category = 20)

gsea %>%
  enrichplot::pairwise_termsim() %>%
  enrichplot::emapplot(cex_label_category = 0.6)

# Biological theme comparison
# --------------------------------------

compareGSEA <- function(ranked_clusters, ...) {
  gsea_list <- lapply(ranked_clusters, function(i) {
    x <- suppressMessages(clusterProfiler::GSEA(i, ...))
    if (class(x) == "gseaResult") {
        as.data.frame(x)
    }
  })

  clusters_levels <- names(ranked_clusters)

  df <- dplyr::bind_rows(gsea_list, .id = "Cluster")
  df[["Cluster"]] <- factor(df$Cluster, levels = clusters_levels)

  new("compareClusterResult",
        compareClusterResult = df,
        geneClusters = ranked_clusters,
        fun = "GSEA",
        .call = match.call(expand.dots = TRUE)
  )
}

#' ep_str_wrap internal string wrapping function
#' @param string the string to be wrapped
#' @param width the maximum number of characters before wrapping to a new line
#' @noRd
ep_str_wrap <- function(string, width) {
    x <- gregexpr(' ', string)
    vapply(seq_along(x),
           FUN = function(i) {
               y <- x[[i]]
               n <- nchar(string[i])
               len <- (c(y,n) - c(0, y)) ## length + 1
               idx <- len > width
               j <- which(!idx)
               if (length(j) && max(j) == length(len)) {
                   j <- j[-length(j)]
               }
               if (length(j)) {
                   idx[j] <- len[j] + len[j+1] > width
               }
               idx <- idx[-length(idx)] ## length - 1
               start <- c(1, y[idx] + 1)
               end <- c(y[idx] - 1, n)
               words <- substring(string[i], start, end)
               paste0(words, collapse="\n")
           },
           FUN.VALUE = character(1)
    )
}

dotplot_for_gsea <- function(
  object,
  show_category = 5,
  font_size = 12,
  label_format = 30,
  title = ""
) {
  df <- as.data.frame(object)

  if (is.null(show_category)) {
    result <- df
  } else if (is.numeric(show_category)) {
    top_n <- function(res, show_category) {
      res %>%
        dplyr::group_split(Cluster) %>%
        purrr::map_dfr(function(df, N) {
          if (length(df$setSize) > N) {
            idx <- order(df$pvalue, decreasing = FALSE)[1:N]
            return(df[idx, ])
          } else {
            return(df)
          }
        }, N = show_category)
    }

    result <- top_n(df, show_category)
  } else {
    result <- subset(df, Description %in% show_category)
  }

  ## remove zero count
  result$Description <- as.character(result$Description) ## un-factor
  GOlevel <- result[, c("ID", "Description")] ## GO ID and Term
  GOlevel <- unique(GOlevel)

  result <- result[result$setSize != 0, ]
  result$Description <- factor(
    result$Description, levels = rev(GOlevel$Description))

  core_count <- function(core_enrichment) {
    core_enrichment %>%
      stringr::str_split("/") %>%
      sapply(length)
  }

  result <- dplyr::mutate(result, Count = core_count(core_enrichment))

  label_func <- function(str) {
    ep_str_wrap(str, label_format)
  }

  if(is.function(label_format)) {
      label_func <- label_format
  }

  p <- ggplot2::ggplot(result, ggplot2::aes_string(
      x = "Cluster", y = "Description", size = "Count")) +
    ggplot2::geom_point(ggplot2::aes_string(color = "p.adjust")) +
    ggplot2::scale_color_continuous(
      low = "red", high = "blue",
      guide = ggplot2::guide_colorbar(reverse=TRUE)) +
    ggplot2::ylab(NULL) +
    ggplot2::ggtitle(title) +
    DOSE::theme_dose(font_size) +
    ggplot2::scale_size_continuous(range = c(3, 8)) +
    ggplot2::scale_y_discrete(labels = label_func)

  p
}

compare_list <- list(
  "VCs_1.1DPI.TR4.MOCK" = prepare_ranked_logFC(
    ds_res$DESeq2, "X1DPI.TR4-X1DPI.MOCK", "Vascular Cells 1"),
  "VCs_1.1DPI.PNR2.MOCK" = prepare_ranked_logFC(
    ds_res$DESeq2, "X1DPI.PNR2-X1DPI.MOCK", "Vascular Cells 1"),
  "VCs_4.1DPI.TR4.MOCK" = prepare_ranked_logFC(
    ds_res$DESeq2, "X1DPI.TR4-X1DPI.MOCK", "Vascular Cells 4"),
  "VCs_4.1DPI.PNR2.MOCK" = prepare_ranked_logFC(
    ds_res$DESeq2, "X1DPI.PNR2-X1DPI.MOCK", "Vascular Cells 4")
)

go_data <- go_data %>%
  dplyr::filter(GO_ID != "")

go2gene <- go_data %>%
  dplyr::select(GO_ID, Gene_ID)
go2name <- go_data %>%
  dplyr::select(GO_ID, GO_Name)

comp_gsea <- compareGSEA(
  compare_list, TERM2GENE = go2gene, TERM2NAME = go2name)

comp_gsea %>%
  dotplot_for_gsea(show_category = 10) +
  ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
dev.off()

# ----------------------------------------

compare_ora <- list(
  "VCs_1.1DPI.TR4.MOCK" = ds_res$DESeq2 %>%
    rhapsodykit::diff_state_significant() %>%
    rhapsodykit::diff_state_pull(
      "X1DPI.TR4-X1DPI.MOCK", "Vascular Cells 1", "gene"),
  "VCs_1.1DPI.PNR2.MOCK" = ds_res$DESeq2 %>%
    rhapsodykit::diff_state_significant() %>%
    rhapsodykit::diff_state_pull(
      "X1DPI.PNR2-X1DPI.MOCK", "Vascular Cells 1", "gene"),
  "VCs_4.1DPI.TR4.MOCK" = ds_res$DESeq2 %>%
    rhapsodykit::diff_state_significant() %>%
    rhapsodykit::diff_state_pull(
      "X1DPI.TR4-X1DPI.MOCK", "Vascular Cells 4", "gene"),
  "VCs_4.1DPI.PNR2.MOCK" = ds_res$DESeq2 %>%
    rhapsodykit::diff_state_significant() %>%
    rhapsodykit::diff_state_pull(
      "X1DPI.PNR2-X1DPI.MOCK", "Vascular Cells 4", "gene")
)

compare_ora <- clusterProfiler::compareCluster(
  compare_ora, fun = "enricher", TERM2GENE = go2gene, TERM2NAME = go2name)

enrichplot::dotplot(compare_ora, showCategory = 10, by = "count") +
  #ggplot2::scale_y_discrete(expand = ggplot2::expansion(c(0, 0), c(1, 1))) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
dev.off()

enrichplot::cnetplot(
  compare_ora, showCategory = 10, node_label = "category", legend_n = 2)
dev.off()
