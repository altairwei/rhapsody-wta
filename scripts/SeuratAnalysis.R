#!/usr/bin/env Rscript

dump_and_quit <- function() {
  message("Dumping frames...")
  dump.frames(to.file = TRUE)
  q(status = 1)
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
  "Seurat",
  "Matrix",
  "readr",
  "ggplot2",
  "cowplot",
  "patchwork",
  "magrittr",
  "gghighlight"
))

if (suppressPackageStartupMessages(!requireNamespace("metap"))) {
  if (!requireNamespace("BiocManager")) install.packages("BiocManager")
  BiocManager::install("multtest")
  BiocManager::install("limma")
  install.packages("metap")
}

#options(error = dump_and_quit)
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
    col_types = readr::cols(
      Cell_Index = readr::col_character()
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


        stim <- unique(unlist(seurat_obj[["stim"]]))
        n_stim <-  length(stim)
        p_umap_split_stim <- Seurat::DimPlot(
          seurat_obj, reduction = "umap", split.by = "stim")
        p_tsne_split_stim <- Seurat::DimPlot(
          seurat_obj, reduction = "tsne", split.by = "stim")
        save_plot(
          file.path(output_folder, "UMAP_Scatter_By_Stim.png"),
          p_umap_split_stim, width = 6 * n_stim
        )
        save_plot(
          file.path(output_folder, "TSNE_Scatter_By_Stim.png"),
          p_tsne_split_stim, width = 6 * n_stim
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

perform_find_all_markers <- function(object, output_folder) {
  markers_df <- Seurat::FindAllMarkers(
    object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  readr::write_csv(markers_df, file.path(output_folder, "Markers_All.csv"))
  top10 <- dplyr::top_n(
    dplyr::group_by(markers_df, cluster), n = 10, wt = avg_logFC)
  marker_heatmap <- Seurat::DoHeatmap(
    object, features = top10$gene) + Seurat::NoLegend()
  save_plot(
    file.path(output_folder, "Markers_Top10_Heatmap.png"),
    marker_heatmap, width = 14, height = 14)
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
    df <- tryCatch(
      Seurat::FindConservedMarkers(
        object, ident.1 = i, grouping.var = "stim"),
      error = function(e) {
        message(toString(e))
        NULL
      }
    )
    if (!is.null(df)) {
      df[["cluster"]] <- i
      df[["gene"]] <- rownames(df)
      rownames(df) <- NULL
      return(df)
    } else {
      return(NULL)
    }
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
      tryCatch(
        {
          diff_genes <- Seurat::FindMarkers(
            ident_cells, ident.1 = couple[2], ident.2 = couple[1])
          diff_genes[["test_group"]] <- couple[2]
          diff_genes[["ctrl_group"]] <- couple[1]
          diff_genes[["comparison"]] <- sprintf(
            "%s_vs_%s", couple[2], couple[1])
          diff_genes[["cluster"]] <- id
          diff_genes[["gene"]] <- rownames(diff_genes)
          diff_genes
        },
        error = function(e) {
          message(toString(e))
          NULL
        }
      )
    })
    do.call(dplyr::bind_rows, diff_list)
  })
  do.call(dplyr::bind_rows, comp_list)
}

#' Plot scatter diagram to compare average expression value of genes on
#'   variable conditions for each cell cluster.
#'
#' @param avg_genes A data frame produced by \code{find_all_avg_expr_genes}
#' @param diff_genes A data frame produced by \code{find_all_diff_expr_genes}
#'
#' @return A named list of ggplot2 object.
plot_avg_expr_genes <- function(avg_genes, diff_genes) {
  diff_genes_by_comp <- split(diff_genes, diff_genes[["comparison"]])
  results <- lapply(diff_genes_by_comp, function(diff_genes_comp) {
    # test condition should be placed on y-axis, ctrl condition on x-axis
    condition_ctrl <- unique(diff_genes_comp[["ctrl_group"]])
    condition_test <- unique(diff_genes_comp[["test_group"]])
    # Calculate coefficient of correlation
    #TODO: fill NA with zero.
    df_by_cluster <- split(avg_genes, avg_genes[["cluster"]])
    xmin <- min(avg_genes[[condition_ctrl]], na.rm = TRUE)
    ymax <- max(avg_genes[[condition_test]], na.rm = TRUE)
    df_by_cluster <- lapply(df_by_cluster, function(x) {
      label_df <- tryCatch(
        {
          pearson_r <- cor(x[[condition_ctrl]], x[[condition_test]],
            method = "pearson", use = "complete.obs")
          data.frame(
            coeff_label = paste(
              "'Pearson ' * italic(R^2) == ", round(pearson_r^2, digits = 2)),
            cluster = unique(x[["cluster"]]),
            position_x = xmin,
            position_y = ymax
          )
        },
        error = function(e) {
          message(
            "Error occurred when processing cluster ",
            unique(x[["cluster"]]), " of ",
            condition_test, "_vs_", condition_ctrl)
          message(toString(e))
          data.frame(
            coeff_label = "",
            cluster = unique(x[["cluster"]]),
            position_x = xmin,
            position_y = ymax
          )
        }
      )
      label_df
    })
    df_coff_label <- do.call(dplyr::bind_rows, df_by_cluster)
    # Highlight DEGs
    diff_genes_to_highlight <- dplyr::filter(diff_genes_comp, p_val < 0.0001)
    df_to_plot <- avg_genes
    df_to_plot$highlight <- "no"
    df_to_plot$highlight[
      df_to_plot$gene %in% diff_genes_to_highlight$gene] <- "yes"
    # Plot facets
    p1 <- ggplot2::ggplot(df_to_plot, ggplot2::aes(
        x = !!ggplot2::sym(condition_ctrl), y = !!ggplot2::sym(condition_test),
        color = highlight, alpha = highlight)) +
      ggplot2::ggtitle(unique(diff_genes_comp[["comparison"]])) +
      ggplot2::geom_point() +
      ggplot2::geom_text(
        data = df_coff_label,
        # Don't inherit color or alpha
        inherit.aes = FALSE,
        mapping = ggplot2::aes(
          x = position_x, y = position_y, label = coeff_label),
        hjust = 0, vjust = 1, size = 8, parse = TRUE) +
      ggplot2::scale_color_manual(values = c("black", "red"), guide = FALSE) +
      ggplot2::scale_alpha_manual(values = c(0.05, 1), guide = FALSE) +
      ggplot2::facet_wrap(ggplot2::vars(cluster)) +
      cowplot::panel_border() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    # Plot whole correlation diagram
    pearson_r_whole <- cor(
      avg_genes[[condition_ctrl]], avg_genes[[condition_test]],
      method = "pearson", use = "complete.obs")
    label_whole <- data.frame(
      coeff_label = paste(
        "'Pearson ' * italic(R^2) == ", round(pearson_r_whole^2, digits = 2)),
      position_x = xmin,
      position_y = ymax
    )
    p2 <- ggplot2::ggplot(avg_genes, ggplot2::aes(
        x = !!ggplot2::sym(condition_ctrl),
        y = !!ggplot2::sym(condition_test))) +
      ggplot2::ggtitle(unique(diff_genes_comp[["comparison"]])) +
      ggplot2::geom_point() +
      ggplot2::geom_text(
        data = label_whole,
        # Don't inherit color or alpha
        inherit.aes = FALSE,
        mapping = ggplot2::aes(
          x = position_x, y = position_y, label = coeff_label),
        hjust = 0, vjust = 1, size = 8, parse = TRUE) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    patchwork::wrap_plots(p1, p2)
  })

  results
}

perform_diff_gene <- function(object, output_folder, draw_plot = TRUE) {
  tryCatch(
    expr = {
      avg_genes_df <- find_all_avg_expr_genes(object)
      readr::write_csv(
        avg_genes_df,
        file.path(output_folder, "DEG_Average_Expression.csv"))
      diff_genes_df <- find_all_diff_expr_genes(object)
      readr::write_csv(
        diff_genes_df, file.path(output_folder, "DEG_All.csv"))
      if (isTRUE(draw_plot)) {
        all_p_list <- plot_avg_expr_genes(avg_genes_df, diff_genes_df)
        for (p_name in names(all_p_list)) {
          save_plot(
            file.path(output_folder,
              sprintf("DEG_Avg_Expr_Scatter_%s.png", p_name)),
            all_p_list[[p_name]],
            width = 32, height = 16
          )
        }
      }
    },
    error = function(e) message(toString(e))
  )
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
