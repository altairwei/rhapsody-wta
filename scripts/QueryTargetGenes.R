#!/usr/bin/env Rscript

library(optparse)
library(magrittr)

parser <- OptionParser()

parser <- add_option(parser,
  c("-g", "--gene-table"),
  dest = "target_genes_file",
  action = "store",
  default = NULL,
  type = "character",
  help = paste0("A CSV file with header, first column must be gene name.")
)

parser <- add_option(parser,
  c("-a", "--annot-col"),
  dest = "annotation_col",
  action = "store",
  default = NULL,
  type = "integer",
  help = paste0("Which column used as annotation.")
)

parser <- add_option(parser,
  c("-c", "--clusters"),
  dest = "selected_clusters",
  action = "store",
  default = NULL,
  type = "character",
  help = paste0("Select cluster to show.")
)

parser <- add_option(parser,
  c("-C", "--output-folder"),
  dest = "output_folder",
  action = "store",
  default = getwd(),
  type = "character",
  help = paste0("Folder to store outputs.")
)

parser <- add_option(parser,
  c("-o", "--csv-table-file"),
  dest = "csvfile",
  action = "store",
  type = "character",
  help = paste0("Write expression table.")
)

parser <- add_option(parser,
  c("-t", "--type"),
  dest = "type",
  action = "store",
  default = "logcounts",
  type = "character",
  help = paste0("Value type of pseudo-bulk data.")
)

parser <- add_option(parser,
  c("-s", "--split-cluster"),
  dest = "split_cluster_into_files",
  action = "store_true",
  default = FALSE,
  help = paste0("Split cluster heatmap into separated files.")
)

parser <- add_option(parser,
  c("--plot"),
  dest = "plot_heatmap",
  action = "store_true",
  default = FALSE,
  help = paste0("Plot heatmap.")
)


arguments <- parse_args2(parser)
options <- arguments$options
options$positionals <- arguments$args

if (length(options$positionals) < 1) {
  stop("At least one position argument is required.\n")
}

stopifnot(
  !is.null(options$target_genes_file)
)

obj <- readRDS(options$positionals[[1]])
genes_table <- readr::read_csv(options$target_genes_file)
pb_arr <- rhapsodykit::make_integrated_pseudo_bulk(obj, options$type)

all_cluster <- dimnames(pb_arr)[[3L]]

if (is.null(options$selected_clusters))
  options$selected_clusters <- all_cluster

if (isTRUE(options$plot_heatmap)) {

  if (!dir.exists(options$output_folder))
    dir.create(options$output_folder)

  if (options$split_cluster_into_files) {
    for (k in all_cluster) {
      png(
        file.path(options$output_folder, sprintf("heatmap_cluster_%s.png", k)),
        units = "px", res = 300,
        width = 400 * ncol(pb_arr),
        height = 70 * nrow(genes_table))

      gene_labels <- genes_table[[1]]
      if (!is.null(options$annotation_col))
        gene_labels <- sprintf("(%s) %s",
          genes_table[[options$annotation_col]], gene_labels)

      p <- rhapsodykit::heatmap_cross_sample(
        pb_arr, genes_table[[1]],
        clusters = k,
        type = options$type,
        show_row_names = TRUE,
        row_labels = gene_labels,
        row_names_side = "left",
        row_dend_side = "right",
        column_names_side = "top",
        column_title = sprintf("Heatmap for Cluster %s", k)
      )

      print(p)

      dev.off()
    }
  } else {
    png(
      file.path(options$output_folder, "heatmap.png"),
      units = "px", res = 300, width = 2000,
      height = 10 * nrow(genes_table) * length(options$selected_clusters), )

    p <- rhapsodykit::heatmap_cross_sample(
      pb_arr, genes_table[[1]],
      clusters = options$selected_clusters,
      type = options$type,
      show_row_names = FALSE,
      row_names_side = "left",
      row_dend_side = "right",
      column_names_side = "top"
    )

    print(p)

    dev.off()
  }
}

if (!is.null(options$csvfile)) {
  pb_small <- pb_arr[genes_table$gene, , , drop = FALSE]
  df_to_write <- purrr::imap(
    purrr::array_branch(pb_small, margin = 3),
    function(mtx, index) {
      df <- data.frame(
        gene = rownames(mtx),
        cluster_id = index,
        mtx,
        row.names = NULL,
        check.names = FALSE
      )

      df
  }) %>% dplyr::bind_rows()

  readr::write_csv(df_to_write, options$csvfile)
}
