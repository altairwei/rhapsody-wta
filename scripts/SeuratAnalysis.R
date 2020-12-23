#!/usr/bin/env Rscript

dump_and_quit <- function() {
  message("Dumping frames...")
  dump.frames(to.file = TRUE)
  q(status = 1)
}

#options(error = dump_and_quit)

if (suppressPackageStartupMessages(!require("Seurat")))
  install.packages("Seurat")
if (suppressPackageStartupMessages(!require("Matrix")))
  install.packages("Matrix")
if (suppressPackageStartupMessages(!require("readr")))
  install.packages("readr")
if (suppressPackageStartupMessages(!require("ggplot2")))
  install.packages("ggplot2")
if (suppressPackageStartupMessages(!require("cowplot")))
  install.packages("cowplot")
if (suppressPackageStartupMessages(!require("patchwork")))
  install.packages("patchwork")
if (suppressPackageStartupMessages(!require("magrittr")))
  install.packages("magrittr")
if (suppressPackageStartupMessages(!require("metap"))) {
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("multtest")
  BiocManager::install("limma")
  install.packages("metap")
}

ggplot2::theme_set(cowplot::theme_cowplot())

save_plot <- function(filename, plot, width = 7, height = 7, ...) {
  message(sprintf("Saving %d x %d in image: %s", width, height, filename))
  ggplot2::ggsave(filename, plot, width = width, height = height, ...)
}

save_png <- function(filename, plot, width = 7, height = 7, dpi = 300) {
  message(sprintf("Saving %d x %d in image: %s", width, height, filename))
  png(filename, width = width, height = height, units = "cm", res = dpi)
  print(plot)
  dev.off()
}

read_raw_csv <- function(base_dir = ".") {
  if (!dir.exists(base_dir)) {
    stop("Directory provided does not exist")
  }

  matrix_loc <- Sys.glob(file.path(base_dir, "*_RSEC_MolsPerCell.csv"))

  if (length(matrix_loc) != 1) {
    stop("`*_RSEC_MolsPerCell.csv` missing or more than one file was found.")
  }

  message(sprintf("Reading %s", matrix_loc))
  df <- readr::read_csv(matrix_loc, comment = "#", progress = TRUE,
    col_types = cols(
      Cell_Index = col_character()
    )
  )

  message("Matrix constructing")
  m <- data.matrix(df[, -1])
  rownames(m) <- df[["Cell_Index"]]
  m <- t(m)
  mtx <- Matrix::Matrix(m)

  rm(df, m)
  gc()

  mtx
}

read_mtx <- function(base_dir = ".") {
  if (!dir.exists(base_dir)) {
    stop("Directory provided does not exist")
  }

  matrix_loc <- Sys.glob(
    file.path(base_dir, "*_Expression_Matrix.mtx"))

  if (length(matrix_loc) != 1) {
    stop("`*_Expression_Matrix.mtx` missing or more than one file was found.")
  }

  colnames_loc <- sprintf("%s.colnames", matrix_loc)
  rownames_loc <- sprintf("%s.rownames", matrix_loc)
  stopifnot(file.exists(colnames_loc, rownames_loc))

  message(sprintf("Reading %s", matrix_loc))
  m <- Matrix::readMM(matrix_loc)
  colnames(m) <- readr::read_lines(colnames_loc)
  rownames(m) <- readr::read_lines(rownames_loc)

  m
}

read_rhapsody_wta <- function(base_dir, use_mtx = FALSE) {
  expr_matrix <- NULL

  if (isTRUE(use_mtx)) {
    expr_matrix <- read_mtx(base_dir)
  } else {
    expr_matrix <- read_raw_csv(base_dir)
  }

  expr_matrix
}

generate_seurat_plots <- function(seurat_obj, output_folder) {
  # Quality Control
  tryCatch(
    expr = {
      p_qc_metrics <- Seurat::VlnPlot(
        seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2
      )
      p_qc_metrics_image_width <- 6

      if (!is.null(seurat_obj@meta.data[["percent.mt"]])) {
        p_qc_metrics <- p_qc_metrics + Seurat::VlnPlot(
          seurat_obj, features = "percent.mt")
        p_qc_metrics_image_width <- p_qc_metrics_image_width + 3
      }

      if (!is.null(seurat_obj@meta.data[["percent.cp"]])) {
        p_qc_metrics <- p_qc_metrics + Seurat::VlnPlot(
          seurat_obj, features = "percent.cp")
        p_qc_metrics_image_width <- p_qc_metrics_image_width + 3
      }
      save_plot(
        file.path(output_folder, "QC_Metrics.png"),
        p_qc_metrics, width = p_qc_metrics_image_width
      )
      save_plot(
        file.path(output_folder, "QC_Count_Feature_Scatter.png"),
        Seurat::FeatureScatter(
          seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      )
    },
    error = function(e) message(toString(e))
  )

  # Feature selection
  tryCatch(
    expr = {
      # Identify the 10 most highly variable genes
      feature_selection_top10 <- head(
        Seurat::VariableFeatures(seurat_obj), 10)
      # plot variable features with and without labels
      p_feature_selection_1 <- Seurat::VariableFeaturePlot(seurat_obj)
      p_feature_selection_2 <- Seurat::LabelPoints(
        plot = p_feature_selection_1,
        points = feature_selection_top10, repel = TRUE)

      save_plot(
        file.path(output_folder, "Feature_Selection.png"),
        p_feature_selection_1 + p_feature_selection_2,
        width = 12
      )
    },
    error = function(e) message(toString(e))
  )

  # PCA Plots
  tryCatch(
    expr = {
      npc <- length(seurat_obj[["pca"]]@stdev)
      p_pca_dim_load <- Seurat::VizDimLoadings(
        seurat_obj, dims = 1:4, reduction = "pca")
      save_plot(
        file.path(output_folder, "PCA_Dim_Loadings.png"),
        p_pca_dim_load
      )
      p_pca_dim_scatter <- Seurat::DimPlot(seurat_obj, reduction = "pca")
      save_plot(
        file.path(output_folder, "PCA_Scatter.png"),
        p_pca_dim_scatter
      )
    },
    error = function(e) message(toString(e))
  )

  tryCatch(
    expr = {
      npc <- length(seurat_obj[["pca"]]@stdev)
      p_jackstraw <- Seurat::JackStrawPlot(seurat_obj, dims = 1:npc) +
        ggplot2::theme(legend.position = "none")
      p_elbow <- Seurat::ElbowPlot(seurat_obj, ndims = npc)
      save_plot(
        file.path(output_folder, "PCA_Dimensionality.png"),
        p_jackstraw + p_elbow,
        width = 12
      )
    },
    error = function(e) message(toString(e))
  )

  tryCatch(
    expr = {
      p_pca_dim_heatmap <- Seurat::DimHeatmap(
        seurat_obj, dims = 1:9, cells = 500, balanced = TRUE)
      save_plot(
        file.path(options$output_folder, "PCA_Dim_Heatmap.png"),
        p_pca_dim_heatmap,
        width = 12
      )
    },
    error = function(e) message(toString(e))
  )

  # UMAP/tSNE Plots
  if (Seurat::DefaultAssay(seurat_obj) == "integrated") {
    tryCatch(
      expr = {
        p_umap_stim <- Seurat::DimPlot(
          seurat_obj, reduction = "umap", group.by = "stim")
        p_umap_cluster <- Seurat::DimPlot(
          seurat_obj, reduction = "umap", label = TRUE)
        save_plot(
          file.path(output_folder, "UMAP_Scatter.png"),
          p_umap_stim + p_umap_cluster,
          width = 12
        )

        p_tsne_stim <- Seurat::DimPlot(
          seurat_obj, reduction = "tsne", group.by = "stim")
        p_tsne_cluster <- Seurat::DimPlot(
          seurat_obj, reduction = "tsne", label = TRUE)
        save_plot(
          file.path(output_folder, "TSNE_Scatter.png"),
          p_tsne_stim + p_tsne_cluster,
          width = 12
        )
      },
      error = function(e) message(toString(e))
    )
  } else {
    tryCatch(
      expr = {
        p_umap_dim <- Seurat::DimPlot(
          seurat_obj, reduction = "umap", label = TRUE)
        save_plot(
          file.path(output_folder, "UMAP_Scatter.png"),
          p_umap_dim
        )
        p_tsne_dim <- Seurat::DimPlot(
          seurat_obj, reduction = "tsne", label = TRUE)
        save_plot(
          file.path(output_folder, "TSNE_Scatter.png"),
          p_tsne_dim
        )
      },
      error = function(e) message(toString(e))
    )
  }

}

single_sample_analysis <- function(
  seurat_obj,
  mt_gene_file = NULL,
  cp_gene_file = NULL,
  dimensionality = 20
) {
  # Quality Control
  if (!is.null(mt_gene_file)) {
    mt_gene_list <- readr::read_lines(mt_gene_file)
    mt_gene_list <- intersect(rownames(seurat_obj[["RNA"]]), mt_gene_list)
    if (length(mt_gene_list) > 0) {
      message(sprintf(
        "Found %d mitochondrial genes in expression matrix.",
        length(mt_gene_list)))
      seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
        seurat_obj, features = mt_gene_list)
    } else {
      message("No mitochondrial genes were found.")
    }
  }

  if (!is.null(cp_gene_file)) {
    cp_gene_list <- readr::read_lines(cp_gene_file)
    cp_gene_list <- intersect(rownames(seurat_obj[["RNA"]]), cp_gene_list)
    if (length(cp_gene_list) > 0) {
      message(sprintf(
        "Found %d chloroplast genes in expression matrix.",
        length(mt_gene_list)))
      seurat_obj[["percent.cp"]] <- Seurat::PercentageFeatureSet(
        seurat_obj, features = cp_gene_list)
    } else {
      message("No chloroplast genes were found.")
    }
  }

  # Normalizing the data
  seurat_obj <- Seurat::NormalizeData(
    seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

  # Feature selection
  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj, selection.method = "vst", nfeatures = 2000)

  # Scaling the data
  all_genes <- rownames(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = all_genes)

  # PCA Plots

  # Perform linear dimensional reduction
  seurat_obj <- Seurat::RunPCA(
    seurat_obj, features = Seurat::VariableFeatures(seurat_obj))
  npc <- length(seurat_obj[["pca"]]@stdev)

  # Determine the 'dimensionality' of the dataset
  seurat_obj <- Seurat::JackStraw(seurat_obj, dims = npc, num.replicate = 100)
  seurat_obj <- Seurat::ScoreJackStraw(seurat_obj, dims = 1:npc)

  # Cluster the cells
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:dimensionality)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)

  # UMAP/tSNE Plots
  #TODO: Make sure umap-learn work properly
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:dimensionality)
  #TODO: Check duplicates manually
  # Workaround: https://github.com/satijalab/seurat/issues/167
  seurat_obj <- Seurat::RunTSNE(
    seurat_obj, dims = 1:dimensionality, check_duplicates = FALSE)

  seurat_obj
}

integrated_sample_analysis <- function(obj_list, dimensionality = 20) {
  stopifnot(is.list(obj_list))

  obj_anchors <- Seurat::FindIntegrationAnchors(
    obj_list, dims = 1:dimensionality)
  obj_combined <- Seurat::IntegrateData(
    anchorset = obj_anchors, dims = 1:dimensionality)
  Seurat::DefaultAssay(obj_combined) <- "integrated"

  # Run the standard workflow for visualization and clustering
  obj_combined <- Seurat::ScaleData(obj_combined)
  obj_combined <- Seurat::RunPCA(obj_combined)
  # t-SNE and Clustering
  #TODO: Make sure umap-learn work properly
  obj_combined <- Seurat::RunUMAP(obj_combined, dims = 1:dimensionality)
  #TODO: Check duplicates manually
  # Workaround: https://github.com/satijalab/seurat/issues/167
  obj_combined <- Seurat::RunTSNE(
    obj_combined, dims = 1:dimensionality, check_duplicates = FALSE)

  obj_combined <- Seurat::FindNeighbors(obj_combined, dims = 1:dimensionality)
  obj_combined <- Seurat::FindClusters(obj_combined, resolution = 0.5)

  obj_combined
}

plot_clustree <- function(object) {
  if (suppressPackageStartupMessages(!require("clustree")))
    install.packages("clustree")
  obj <- Seurat::FindClusters(
    object, resolution = seq(0.4, 1.6, 0.2))
  clustree::clustree(obj@meta.data, prefix = "RNA_snn_res.")
}

find_all_conserved_markers <- function(object) {
  Seurat::DefaultAssay(object) <- "RNA"
  idents_all <- sort(unique(Seurat::Idents(object)))
  df_list <- lapply(idents_all, function(i) {
    df <- Seurat::FindConservedMarkers(
      object, ident.1 = i, grouping.var = "stim")
    df[["cluster"]] <- i
    df[["gene"]] <- rownames(df)
    rownames(df) <- NULL
    df
  })
  do.call(dplyr::bind_rows, df_list)
}

find_all_avg_expr_genes <- function(object) {
  idents_all <- sort(unique(Seurat::Idents(object)))
  df_list <- lapply(idents_all, function(i) {
    ident_cells <- subset(object, idents = i)
    # Specify identity of cells based on value of meta.data[["stim"]]
    Seurat::Idents(ident_cells) <- "stim"
    # AverageExpression will be applied to every `Assay` object.
    avg_ident_cells <- log1p(Seurat::AverageExpression(ident_cells)$RNA)
    avg_ident_cells[["cluster"]] <- i
    avg_ident_cells[["gene"]] <- rownames(avg_ident_cells)
    avg_ident_cells
  })
  do.call(dplyr::bind_rows, df_list)
}

#' Find differential expressed genes on variable conditions
#'  for each cell cluster.
#'
#' @param object Seurat object.
#'
#' @return A data frame. 
find_all_diff_expr_genes <- function(object) {
  idents_all <- sort(unique(Seurat::Idents(object)))
  stim <- unique(unlist(object[["stim"]]))
  # Produce an permutation of all conditions
  all_comb <- as.data.frame(combn(stim, 2), stringsAsFactors = FALSE)
  names(all_comb) <- NULL
  # Find diff genes on each comparision
  comp_list <- lapply(all_comb, function(couple) {
    message(sprintf("Calculating comparision %s vs %s", couple[2], couple[1]))
    # Perform on each cell cluster
    diff_list <- lapply(idents_all, function(id) {
      ident_cells <- subset(object, idents = id)
      # Specify identity of cells based on value of meta.data[["stim"]]
      Seurat::Idents(ident_cells) <- "stim"
      message("Calculating cluster ", id)
      # We assume that the control group is the first one
      #   and the stimulated group is the second one.
      diff_genes <- Seurat::FindMarkers(
        ident_cells, ident.1 = couple[2], ident.2 = couple[1])
      diff_genes[["comparison"]] <- sprintf("%s_vs_%s", couple[2], couple[1])
      diff_genes[["cluster"]] <- id
      diff_genes[["gene"]] <- rownames(diff_genes)
      diff_genes
    })
    do.call(dplyr::bind_rows, diff_list)
  })
  do.call(dplyr::bind_rows, comp_list)
}

#' Plot scatter diagram to compare average expression value of genes on
#'   variable conditions for each cell cluster.
#'
#' @param df A data frame produced by \code{find_all_avg_expr_genes}
#'
#' @return A list of ggplot2 object.
plot_avg_expr_genes <- function(df) {
  stim <- names(df)[!names(df) %in% c("cluster", "gene")]
  # Produce an permutation of all conditions
  all_comb <- as.data.frame(combn(stim, 2), stringsAsFactors = FALSE)
  names(all_comb) <- NULL
  lapply(all_comb, function(couple) {
    ggplot2::ggplot(df, ggplot2::aes(
        !!ggplot2::sym(couple[1]), !!ggplot2::sym(couple[2]))) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(vars(cluster)) +
      cowplot::panel_border()
  })
}

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
    integrate = FALSE
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
    markers_df <- Seurat::FindAllMarkers(
      obj_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    readr::write_csv(markers_df, file.path(output_folder, "Markers_All.csv"))

    # Conserved markers across conditions
    markers_conserved_df <- find_all_conserved_markers(obj_combined)
    readr::write_csv(
      markers_conserved_df, file.path(output_folder, "Markers_Conserved.csv"))

    # Identify differential expressed genes across conditions
    tryCatch(
      expr = {
        avg_genes_df <- find_all_avg_expr_genes(obj_combined)
        readr::write_csv(
          avg_genes_df,
          file.path(output_folder, "DEG_Average_Expression.csv"))
        diff_genes_df <- find_all_diff_expr_genes(obj_combined)
        readr::write_csv(
          diff_genes_df, file.path(output_folder, "DEG_All.csv"))
        if (isTRUE(options$draw_plot)) {
          save_plot(
            file.path(output_folder, "DEG_Avg_Expr_Scatter.png"),
            plot_avg_expr_genes(avg_genes_df),
            width = 16, height = 16
          )
        }
      },
      error = function(e) message(toString(e))
    )

    if (isTRUE(options$produce_cache)) {
      saveRDS(obj_combined,
        file = file.path(output_folder, "Seurat_Object_Combined.rds"))
    }

    if (isTRUE(options$draw_plot)) {
      generate_seurat_plots(obj_combined, output_folder)
    }
  } else {
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

      markers_df <- Seurat::FindAllMarkers(
        seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      readr::write_csv(markers_df, file.path(output_folder, "Markers_All.csv"))

      if (isTRUE(options$produce_cache)) {
        saveRDS(seurat_obj,
          file = file.path(output_folder, "Seurat_Object.rds"))
      }

      if (isTRUE(options$draw_plot)) {
        generate_seurat_plots(seurat_obj, output_folder)
      }

    })
  }

}
