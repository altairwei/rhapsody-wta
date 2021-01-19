#!/usr/bin/env Rscript

load_sources <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  script_arg_prefix <- "--file="
  script_name <- sub(
    script_arg_prefix, "", full_args[grep(script_arg_prefix, full_args)])
  real_name <- Sys.readlink(script_name)
  if (isTRUE(nzchar(real_name, keepNA = TRUE))) {
    script_name <- real_name
  }
  script_dir <- dirname(script_name)
  deps_dir <- file.path(script_dir, "RAnalysis")
  devtools::load_all(deps_dir)
}

require_dependencies <- function(x) {
  stopifnot(is.character(x))
  dep_status <- suppressPackageStartupMessages(
    sapply(x, requireNamespace)
  )
  pkg_to_install <- names(dep_status[dep_status == FALSE])
  install.packages(pkg_to_install)
}

require_dependencies(c(
  "optparse",
  "devtools"
))

load_sources()

ggplot2::theme_set(cowplot::theme_cowplot())

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  options <- list(
    positionals = character(0),
    draw_plot = FALSE,
    compress = FALSE,
    use_mtx = FALSE,
    use_cache = NULL,
    produce_cache = FALSE,
    output_folder = getwd(),
    mt_gene_file = NULL,
    cp_gene_file = NULL,
    integrate = FALSE,
    process = 1L
  )

  optind <- 1
  while (optind <= length(args)) {
    switch(args[optind],
      "--plot" = {
        options$draw_plot <- TRUE
      },
      "-c" = {
        # Convert expression csv file produced by Rhapsody WTA pipeline
        #   into MatrixMarket format.
        options$compress <- TRUE
      },
      "-m" = {
        # Try to find MatrixMarket format instead of csv format table.
        options$use_mtx <- TRUE
      },
      "-O" = {
        optind <- optind + 1
        options$output_folder <- args[optind]
      },
      "--mt-list" = {
        # Mitochondrial gene list
        optind <- optind + 1
        options$mt_gene_file <- args[optind]
      },
      "--cp-list" = {
        # Chloroplast gene list
        optind <- optind + 1
        options$cp_gene_file <- args[optind]
      },
      "--use-cache" = {
        optind <- optind + 1
        options$use_cache <- args[optind]
      },
      "--produce-cache" = {
        options$produce_cache <- TRUE
      },
      "--integrate" = {
        options$integrate <- TRUE
      },
      "--process" = {
        optind <- optind + 1
        options$process <- as.integer(args[optind])
      },
      {
        if (startsWith(args[optind], "-")) {
          stop(sprintf("Unknown option: %s", args[optind]))
        } else {
          options$positionals <- append(options$positionals, args[optind])
        }
      }
    )
    optind <- optind + 1
  }

  if (length(options$positionals) < 1) {
    stop("At least one position argument is required.\n")
  }

  if (options$process > 1) {
    if (suppressPackageStartupMessages(!requireNamespace("future")))
      install.packages("future")
    future::plan("multiprocess", workers = options$process)
    # Set global size to 2GB
    options(future.globals.maxSize = 4 * 1024^3)
  }

  if (isTRUE(options$integrate)) {
    obj_combined <- NULL
    if (!is.null(options$use_cache)) {
      obj_combined <- readRDS(options$use_cache)
    } else {
      obj_list <- lapply(options$positionals, function(base_dir) {
        expr_matrix <- read_rhapsody_wta(base_dir, options$use_mtx)
        seurat_obj <- Seurat::CreateSeuratObject(
          counts = expr_matrix, project = basename(base_dir))
        seurat_obj$stim <- basename(base_dir)
        # TODO: 增加 SCTransform 的选项！
        seurat_obj <- Seurat::NormalizeData(seurat_obj)
        seurat_obj <- Seurat::FindVariableFeatures(
          seurat_obj, selection.method = "vst", nfeatures = 2000)
        seurat_obj
      })

      obj_combined <- integrated_sample_analysis(obj_list)
    }

    output_folder <- options$output_folder
    if (!dir.exists(output_folder)) {
      dir.create(output_folder)
    }

    #FIXME: Does this make sense?
    perform_find_all_markers(obj_combined, output_folder)

    # Conserved markers across conditions
    markers_conserved_df <- find_all_conserved_markers(obj_combined)
    readr::write_csv(
      markers_conserved_df, file.path(output_folder, "Markers_Conserved.csv"))

    # Identify differential expressed genes across conditions
    perform_diff_gene(obj_combined, output_folder, options$draw_plot)

    if (isTRUE(options$produce_cache)) {
      saveRDS(obj_combined,
        file = file.path(output_folder, "Seurat_Object_Combined.rds"))
    }

    if (isTRUE(options$draw_plot)) {
      generate_seurat_plots(obj_combined, output_folder)
    }
  } else {
    invisible(
      lapply(options$positionals, function(base_dir) {
        expr_matrix <- read_rhapsody_wta(base_dir, options$use_mtx)

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

        seurat_obj <- single_sample_analysis(
          seurat_obj, options$mt_gene_file, options$cp_gene_file)

        perform_find_all_markers(seurat_obj, output_folder)

        if (isTRUE(options$produce_cache)) {
          saveRDS(seurat_obj,
            file = file.path(output_folder, "Seurat_Object.rds"))
        }

        if (isTRUE(options$draw_plot)) {
          generate_seurat_plots(seurat_obj, output_folder)
        }

      })
    )
  }

}
