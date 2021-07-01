library(magrittr)

contrast <- ""
cluster <- ""

# Over Representation Analysis
# --------------------------------------

ora <- rhapsodykit::enrich_perform_ora(
  ds_res$DESeq2, go_data,
  contrast,
  cluster
)

fc <- rhapsodykit::enrich_significant_lfc(
  ds_res$DESeq2, contrast, cluster)

ora %>%
  rhapsodykit::enrich_barplot()
dev.off()

ora %>%
  rhapsodykit::enrich_cnetplot(fold_change = fc, show_category = 20)
dev.off()

ora %>%
  enrichplot::pairwise_termsim() %>%
  enrichplot::emapplot(cex_label_category = 0.6)
dev.off()

# Gene Set Enrichment Analysis
# --------------------------------------

gsea <- rhapsodykit::enrich_perform_gsea(
  ds_res$DESeq2, go_data,
  contrast, cluster
)

gsea %>%
  rhapsodykit::enrich_cnetplot(
    fold_change = gsea@geneList, show_category = 20)
dev.off()

gsea %>%
  enrichplot::pairwise_termsim() %>%
  enrichplot::emapplot(cex_label_category = 0.6)
dev.off()

# Biological theme comparison for ORA
# ----------------------------------------

go_data <- go_data %>%
  dplyr::filter(GO_ID != "")

go2gene <- go_data %>%
  dplyr::select(GO_ID, Gene_ID)
go2name <- go_data %>%
  dplyr::select(GO_ID, GO_Name)

prepare_data_list <- function(
  data,
  contrasts,
  clusters,
  fun
) {
  res <- purrr::imap(clusters, function(clr, clr_name) {
    purrr::imap(contrasts, function(con, con_name) {
      fun(data, con, clr)
    })
  })

  unlist(res, recursive = FALSE)
}

compare_ora <- prepare_data_list(
  ds_res$DESeq2,
  contrasts = c(
    "3DPI.TR4.MOCK" = "X3DPI.TR4-X3DPI.MOCK",
    "3DPI.PNR2.MOCK" = "X3DPI.PNR2-X3DPI.MOCK"
  ),
  clusters = c(
    "VCs_6" = "Vascular Cells 6",
    "VCs_5" = "Vascular Cells 5",
    "VCs_1" = "Vascular Cells 1",
    "VCs_2" = "Vascular Cells 2",
    "VCs_4" = "Vascular Cells 4",
    "ECs" = "Epidermal Cells"
  ),
  fun = function(data, con, clr) {
    rhapsodykit::enrich_significant_lfc(
      data, con, clr, 0.01, 2) %>% names()
  }
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

compare_ora %>%
  enrichplot::pairwise_termsim() %>%
  enrichplot::emapplot(legend_n = 2)
dev.off()

# Biological theme comparison for GSEA
# --------------------------------------

compare_list <- prepare_data_list(
  ds_res$DESeq2,
  contrasts = c(
    "3DPI.TR4.MOCK" = "X3DPI.TR4-X3DPI.MOCK",
    "3DPI.PNR2.MOCK" = "X3DPI.PNR2-X3DPI.MOCK"
  ),
  clusters = c(
    "VCs_6" = "Vascular Cells 6",
    "VCs_5" = "Vascular Cells 5",
    "VCs_1" = "Vascular Cells 1",
    "VCs_2" = "Vascular Cells 2",
    "VCs_4" = "Vascular Cells 4",
    "ECs" = "Epidermal Cells"
  ),
  fun = function(data, con, clr) {
    rhapsodykit::enrich_ranked_lfc(data, con, clr)
  }
)

comp_gsea <- rhapsodykit::enrich_compare_gsea(
  compare_list, TERM2GENE = go2gene, TERM2NAME = go2name)

comp_gsea %>%
  rhapsodykit::dotplot_for_gsea(show_category = 10) +
  ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
dev.off()
