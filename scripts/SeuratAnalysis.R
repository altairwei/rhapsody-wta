#!/usr/bin/env Rscript

library(optparse)

ggplot2::theme_set(cowplot::theme_cowplot())

parser <- OptionParser()
parser <- add_option(parser,
  c("--analysis"),
  dest = "analysis",
  action = "store",
  default = "standalone",
  help = paste0("Which analysis should be applied to samples.",
    "[default: %default]"))
parser <- add_option(parser,
  c("-i", "--integrate"),
  dest = "analysis",
  action = "callback",
  callback = function(...) "integrate",
  help = paste0("Integrate data by using Seurat."))
parser <- add_option(parser,
  c("--merge"),
  dest = "analysis",
  action = "callback",
  callback = function(...) "merge",
  help = paste0("Merge all samples into one object without batch effect ",
    "correction."))
parser <- add_option(parser,
  c("-d", "--draw-plot"),
  dest = "draw_plot",
  action = "store_true",
  default = FALSE,
  type = "logical",
  help = paste0("Draw plots for Seurat analysis. [default: %default]"))
parser <- add_option(parser,
  c("--diff-gene"),
  dest = "diff_gene",
  action = "store_true",
  default = FALSE,
  type = "logical",
  help = paste0("Perform differential expression analysis to ",
    "find markers. [default: %default]"))
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
  help = paste0("How many neighbors to use when picking anchors.",
    " [default: %default]")
)
parser <- add_option(parser,
  c("-r", "--reduction"),
  dest = "reduction",
  action = "store",
  default = "cca",
  type = "character",
  help = paste0("Dimensional reduction to perform when finding anchors.",
    " [default: %default]")
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
parser <- add_option(parser,
  c("-L", "--data-list"),
  dest = "data_list",
  action = "store",
  default = NULL,
  type = "character",
  help = paste0("Provide a CSV file which lists data folder in `data_folder`",
    "column and reference in `use_as_ref` column."))

arguments <- parse_args2(parser)
options <- arguments$options
options$positionals <- arguments$args

if (is.character(options$reference)) {
  options$reference <- type.convert(
    strsplit(options$reference, ",")[[1]], as.is = TRUE)
}

data_folders <- character(0)

if (!is.null(options$data_list)) {
  data_df <- readr::read_csv(
    options$data_list,
    col_names = TRUE,
    col_types = readr::cols(
      data_folder = readr::col_character(),
      use_as_ref = readr::col_logical()
    )
  )

  options$data_folders <- ifelse(
    startsWith(data_df$data_folder, "/"),
    data_df$data_folder,
    file.path(dirname(options$data_list), data_df$data_folder)
  )

  if (is.null(options$reference) &&
      length(which(data_df$use_as_ref)) != 0) {
    options$reference <- which(data_df$use_as_ref)
  }

} else {
  options$data_folders <- options$positionals
}

if (length(options$data_folders) < 1) {
  stop("At least one position argument is required.\n")
}

message("Options:")
message(
  paste(
    capture.output(
      str(options, vec.len = 50, no.list = TRUE)
    ),
    collapse = "\n"
  )
)

if (options$process > 1) {
  if (suppressPackageStartupMessages(!requireNamespace("future")))
    install.packages("future")
  future::plan("multiprocess", workers = options$process)
  options(future.globals.maxSize = options$memory * 1024^3)
}

run_analysis <- switch(options$analysis,

  integrate = function() {
    #------------------------------------------------------------
    # Integrate all sample into one object.
    #------------------------------------------------------------

    obj_combined <- NULL
    if (!is.null(options$use_cache)) {
      obj_combined <- readRDS(options$use_cache)
    } else {
      # Convert raw data into Seurat object for each sample.
      obj_list <- lapply(options$data_folders, function(base_dir) {
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

      # Integrate all samples.
      obj_combined <- rhapsodykit::integrated_sample_analysis(
        obj_list, reduction = options$reduction, k.anchor = options$anchors,
        reference = options$reference
      )
    }

    output_folder <- options$output_folder
    if (!dir.exists(output_folder)) {
      dir.create(output_folder)
    }

    if (isTRUE(options$diff_gene)) {
      #FIXME: Does this make sense?
      rhapsodykit::perform_find_all_markers(obj_combined, output_folder)

      # Conserved markers across conditions
      markers_conserved_df <- rhapsodykit::find_all_conserved_markers(
        obj_combined)
      readr::write_csv(
        markers_conserved_df, file.path(output_folder, "Markers_Conserved.csv"))

      # Identify differential expressed genes across conditions
      rhapsodykit::perform_diff_gene(
        obj_combined, output_folder, options$draw_plot)
    }

    if (isTRUE(options$produce_cache)) {
      saveRDS(obj_combined,
        file = file.path(output_folder, "Seurat_Object_Combined.rds"))
    }

    if (isTRUE(options$draw_plot)) {
      rhapsodykit::generate_seurat_plots(obj_combined, output_folder)
    }
  },

  merge = function() {
    # Convert raw data into Seurat object for each sample.
    obj_list <- lapply(options$data_folders, function(base_dir) {
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

    output_folder <- options$output_folder
    if (!dir.exists(output_folder)) {
      dir.create(output_folder)
    }

    obj_merged <- rhapsodykit::merge_sample(obj_list)

    if (isTRUE(options$produce_cache)) {
      saveRDS(obj_merged,
        file = file.path(output_folder, "Seurat_Object_Merged.rds"))
    }
  },

  standalone = function() {
    #------------------------------------------------------------
    # Process all samples in parallel.
    #------------------------------------------------------------

    invisible(
      lapply(options$data_folders, function(base_dir) {
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

        if (isTRUE(options$diff_gene)) {
          rhapsodykit::perform_find_all_markers(seurat_obj, output_folder)
        }

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
)

run_analysis()