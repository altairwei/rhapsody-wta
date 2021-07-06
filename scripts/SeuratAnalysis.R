#!/usr/bin/env Rscript

require_dependencies <- function(x) {
  stopifnot(is.character(x))
  dep_status <- suppressPackageStartupMessages(
    sapply(x, requireNamespace)
  )
  pkg_to_install <- names(dep_status[dep_status == FALSE])
  if (length(pkg_to_install) > 0)
    install.packages(pkg_to_install)
}

require_dependencies(c(
  "optparse",
  "rhapsodykit"
))

library(optparse)

ggplot2::theme_set(cowplot::theme_cowplot())

parser <- OptionParser()
parser <- add_option(parser,
  c("-d", "--draw-plot"),
  dest = "draw_plot",
  action = "store_true",
  default = FALSE,
  type = "logical",
  help = paste0("Draw plots for Seurat analysis. [default: %default]"))
parser <- add_option(parser,
  c("--compress"),
  dest = "compress",
  action = "store_true",
  default = FALSE,
  type = "logical",
  help = paste0("Convert expression csv file produced by Rhapsody WTA",
    " pipeline into MatrixMarket format. [default: %default]"))
parser <- add_option(parser,
  c("-m", "--use-matrix"),
  dest = "use_mtx",
  action = "store_true",
  default = FALSE,
  type = "logical",
  help = paste0("Try to find MatrixMarket format instead of csv format",
    " table. [default: %default]"))
parser <- add_option(parser,
  c("-c", "--produce-cache"),
  dest = "produce_cache",
  action = "store_true",
  default = FALSE,
  type = "logical",
  help = paste0("Produce Seurat obj by `saveRDS`. [default: %default]"))
parser <- add_option(parser,
  c("--use-cache"),
  dest = "use_cache",
  action = "store",
  default = NULL,
  type = "character",
  help = paste0("Use cache file produced by `saveRDS`. [default: %default]"))
parser <- add_option(parser,
  c("-O", "--output-folder"),
  dest = "output_folder",
  action = "store",
  default = getwd(),
  type = "character",
  help = paste0("Folder to save outputs. [default: current working dir]"))
parser <- add_option(parser,
  c("--cp-list"),
  dest = "cp_gene_file",
  action = "store",
  default = NULL,
  type = "character",
  help = paste0("Chloroplast gene list. [default: %default]"))
parser <- add_option(parser,
  c("--mt-list"),
  dest = "mt_gene_file",
  action = "store",
  default = NULL,
  type = "character",
  help = paste0("Mitochondrial gene list. [default: %default]"))
parser <- add_option(parser,
  c("-i", "--integrate"),
  dest = "integrate",
  action = "store_true",
  default = FALSE,
  type = "logical",
  help = paste0("Integrate data by useing Seurat. [default: %default]"))
parser <- add_option(parser,
  c("--group-rule"),
  dest = "group_rule",
  action = "store",
  default = NULL,
  type = "character",
  help = paste0("The pattern to extract group information from sample name.")
)
parser <- add_option(parser,
  c("-R", "--reference"),
  dest = "reference",
  action = "store",
  default = NULL,
  type = "character",
  help = paste0("Perform reference-based integration based on provided",
    "reference locations (integer). Multiple references should be ",
    "separated by comma.")
)
parser <- add_option(parser,
  c("--anchors"),
  dest = "anchors",
  action = "store",
  default = 5,
  type = "integer",
  help = paste0("How many neighbors to use when picking anchors. [default: %default]")
)
parser <- add_option(parser,
  c("-r", "--reduction"),
  dest = "reduction",
  action = "store",
  default = "cca",
  type = "character",
  help = paste0("Dimensional reduction to perform when finding anchors. [default: %default]")
)
parser <- add_option(parser,
  c("--process"),
  dest = "process",
  action = "store",
  default = 1L,
  type = "integer",
  help = paste0("How many processes to use. [default: %default]"))
parser <- add_option(parser,
  c("--mem-size"),
  dest = "memory",
  action = "store",
  default = 4L,
  type = "integer",
  help = paste0("How many memory (GB) to use. [default: %default GB]"))
arguments <- parse_args2(parser)
options <- arguments$options
options$positionals <- arguments$args

if (is.character(options$reference)) {
  options$reference <- type.convert(
    strsplit(options$reference, ",")[[1]], as.is = TRUE)
}

if (length(options$positionals) < 1) {
  stop("At least one position argument is required.\n")
}

if (options$process > 1) {
  if (suppressPackageStartupMessages(!requireNamespace("future")))
    install.packages("future")
  future::plan("multiprocess", workers = options$process)
  options(future.globals.maxSize = options$memory * 1024^3)
}

if (isTRUE(options$integrate)) {
  obj_combined <- NULL
  if (!is.null(options$use_cache)) {
    obj_combined <- readRDS(options$use_cache)
  } else {
    obj_list <- lapply(options$positionals, function(base_dir) {
      expr_matrix <- rhapsodykit::read_rhapsody_wta(base_dir, options$use_mtx)
      seurat_obj <- Seurat::CreateSeuratObject(
        counts = expr_matrix, project = basename(base_dir))
      seurat_obj$sample <- basename(base_dir)

      if (!is.null(options$group_rule)) {
        seurat_obj$group <- stringr::str_extract(
          basename(base_dir), options$group_rule)
      } else {
        seurat_obj$group <- seurat_obj$sample
      }

      seurat_obj
    })

    obj_combined <- rhapsodykit::integrated_sample_analysis(
      obj_list, reduction = options$reduction, k.anchor = options$anchors,
      reference = options$reference
    )
  }

  output_folder <- options$output_folder
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }

  #FIXME: Does this make sense?
  rhapsodykit::perform_find_all_markers(obj_combined, output_folder)

  # Conserved markers across conditions
  markers_conserved_df <- rhapsodykit::find_all_conserved_markers(obj_combined)
  readr::write_csv(
    markers_conserved_df, file.path(output_folder, "Markers_Conserved.csv"))

  # Identify differential expressed genes across conditions
  rhapsodykit::perform_diff_gene(obj_combined, output_folder, options$draw_plot)

  if (isTRUE(options$produce_cache)) {
    saveRDS(obj_combined,
      file = file.path(output_folder, "Seurat_Object_Combined.rds"))
  }

  if (isTRUE(options$draw_plot)) {
    rhapsodykit::generate_seurat_plots(obj_combined, output_folder)
  }
} else {
  invisible(
    lapply(options$positionals, function(base_dir) {
      expr_matrix <- rhapsodykit::read_rhapsody_wta(base_dir, options$use_mtx)

      if (isTRUE(options$compress) && !isTRUE(options$use_mtx)) {
        input_file <- Sys.glob(
          file.path(base_dir, "*_RSEC_MolsPerCell.csv"))
        mtx_file <- sub(
          "_RSEC_MolsPerCell.csv$",
          "_Expression_Matrix.mtx", input_file)
        colnames_file <- sprintf("%s.colnames", mtx_file)
        rownames_file <- sprintf("%s.rownames", mtx_file)

        message(sprintf("Writting %s", mtx_file))
        Matrix::writeMM(expr_matrix, file = mtx_file)
        write(colnames(expr_matrix), colnames_file)
        write(rownames(expr_matrix), rownames_file)
      }

      output_folder <- file.path(base_dir, "SeuratAnalysis")

      if (!dir.exists(output_folder)) {
        dir.create(output_folder)
      }

      seurat_obj <- Seurat::CreateSeuratObject(
        counts = expr_matrix, project = basename(base_dir))

      seurat_obj <- rhapsodykit::single_sample_analysis(
        seurat_obj, options$mt_gene_file, options$cp_gene_file)

      rhapsodykit::perform_find_all_markers(seurat_obj, output_folder)

      if (isTRUE(options$produce_cache)) {
        saveRDS(seurat_obj,
          file = file.path(output_folder, "Seurat_Object.rds"))
      }

      if (isTRUE(options$draw_plot)) {
        rhapsodykit::generate_seurat_plots(seurat_obj, output_folder)
      }

    })
  )
}
