#!/usr/bin/env Rscript

dump_and_quit <- function() {
  message("Dumping frames...")
  dump.frames(to.file = TRUE)
  q(status = 1)
}

#options(error = dump_and_quit)

library(Seurat)
library(Matrix)
library(readr)
library(ggplot2)
library(patchwork)

# TODO: 将数据处理过程与画图分开，然后缓冲处理成功了数据对象。

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

read_rhapsody_wta <- function(base_dir = ".") {
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

generate_seurat_plots <- function(seurat_obj, output_dir) {
  # Quality Control
  p_qc_metrics <- Seurat::VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA"),
    ncol = 2
  )
  p_qc_metrics_image_width <- 6

  if (!is.null(seurat_obj[["percent.mt"]])) {
    p_qc_metrics <- p_qc_metrics + Seurat::VlnPlot(
      seurat_obj, features = "percent.mt")
    p_qc_metrics_image_width <- p_qc_metrics_image_width + 3
  }

  if (!is.null(seurat_obj[["percent.cp"]])) {
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

  # Feature selection
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

  # PCA Plots
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
  # p_pca_dim_heatmap <- Seurat::DimHeatmap(
  #   seurat_obj, dims = 1, cells = 500, balanced = TRUE)
  # save_plot(
  #   file.path(options$output_folder, "PCA_Dim_Heatmap.png"),
  #   p_pca_dim_heatmap,
  #   width = 12
  # )
  p_jackstraw <- Seurat::JackStrawPlot(seurat_obj, dims = 1:npc) +
    ggplot2::theme(legend.position = "none")
  p_elbow <- Seurat::ElbowPlot(seurat_obj, ndims = npc)
  save_plot(
    file.path(output_folder, "PCA_Dimensionality.png"),
    p_jackstraw + p_elbow,
    width = 12
  )

  # UMAP/tSNE Plots
  p_umap_dim <- Seurat::DimPlot(seurat_obj, reduction = "umap")
  save_plot(
    file.path(output_folder, "UMAP_Scatter.png"),
    p_umap_dim
  )
  p_tsne_dim <- Seurat::DimPlot(seurat_obj, reduction = "tsne")
  save_plot(
    file.path(output_folder, "TSNE_Scatter.png"),
    p_tsne_dim
  )
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  options <- list(
    positionals = character(0),
    draw_plot = FALSE,
    compress = FALSE,
    use_mtx = FALSE,
    use_cache_file = NULL,
    cache = FALSE,
    output_folder = NULL,
    mt_gene_file = NULL,
    cp_gene_file = NULL
  )

  optind <- 1
  while (optind <= length(args)) {
    switch(args[optind],
      "-p" = {
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
      "--use-cache-file" = {
        optind <- optind + 1
        options$use_cache_file <- args[optind]
      },
      "--cache" = {
        options$cache <- TRUE
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

  options$base_dir <- options$positionals[1]

  # Reading expression matrix
  expr_matrix <- NULL
  if (isTRUE(options$use_mtx)) {
    expr_matrix <- read_mtx(options$base_dir)
  } else {
    expr_matrix <- read_rhapsody_wta(options$base_dir)
  }

  if (isTRUE(options$compress) && !isTRUE(options$use_mtx)) {
    input_file <- Sys.glob(
      file.path(options$base_dir, "*_RSEC_MolsPerCell.csv"))
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

  if (is.null(options$output_folder)) {
    options$output_folder <- file.path(options$base_dir, "SeuratAnalysis")
  }

  if (!dir.exists(options$output_folder)) {
    dir.create(options$output_folder)
  }

  seurat_obj <- Seurat::CreateSeuratObject(
    counts = expr_matrix, project = basename(options$base_dir))

  #########################
  # Quality Control
  #########################

  if (!is.null(options$mt_gene_file)) {
    mt_gene_list <- readr::read_lines(options$mt_gene_file)
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

  if (!is.null(options$cp_gene_file)) {
    cp_gene_list <- readr::read_lines(options$cp_gene_file)
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

  #########################
  # Normalizing the data
  #########################

  seurat_obj <- Seurat::NormalizeData(
    seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

  #########################
  # Feature selection
  #########################

  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj, selection.method = "vst", nfeatures = 2000)

  #########################
  # Scaling the data
  #########################

  all_genes <- rownames(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = all_genes)

  #########################
  # PCA Plots
  #########################

  # Perform linear dimensional reduction
  seurat_obj <- Seurat::RunPCA(
    seurat_obj, features = Seurat::VariableFeatures(seurat_obj))
  npc <- length(seurat_obj[["pca"]]@stdev)

  # Determine the 'dimensionality' of the dataset
  seurat_obj <- Seurat::JackStraw(seurat_obj, dims = npc, num.replicate = 100)
  seurat_obj <- Seurat::ScoreJackStraw(seurat_obj, dims = 1:npc)

  #########################
  # Cluster the cells
  #########################

  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:10)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)

  #########################
  # UMAP/tSNE Plots
  #########################

  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:10)
  seurat_obj <- Seurat::RunTSNE(seurat_obj, dims = 1:10)

  if (isTRUE(options$cache)) {
    saveRDS(pbmc, file = file.path(
      options$base_dir, "SeuratAnalysis", "Seurat_Object.rds"))
  }
}
