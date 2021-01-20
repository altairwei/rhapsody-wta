#!/usr/bin/env Rscript

require_dependencies <- function(x) {
  stopifnot(is.character(x))
  dep_status <- suppressPackageStartupMessages(
    sapply(x, requireNamespace)
  )
  pkg_to_install <- names(dep_status[dep_status == FALSE])
  install.packages(pkg_to_install)
}

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

require_dependencies(c(
  "optparse",
  "devtools",
  "magrittr"
))

library(magrittr)
library(optparse)

load_sources()

options(stringsAsFactors = FALSE)

parser <- OptionParser()
parser <- add_option(parser,
  c("-m", "--method"),
  dest = "method",
  action = "store",
  default = "avg",
  type = "character",
  help = paste0("Create pseudo-bulk data by summing/averaging",
    " the counts from cells. Choices are `avg` and `sum`. [default: %default]"))
parser <- add_option(parser,
  c("--sum"),
  dest = "method",
  action = "callback",
  callback = function(...) "sum",
  help = "Create pseudo-bulk data by summing the counts from cells.")
parser <- add_option(parser,
  c("--avg"),
  dest = "method",
  action = "callback",
  callback = function(...) "avg",
  help = "Create pseudo-bulk data by averaging the counts from cells.")
parser <- add_option(parser,
  c("--sctransform"),
  action = "store_true",
  default = FALSE,
  help = "Use Seurat::SCTransform to perform normalization.")
parser <- add_option(parser,
  c("-O", "--output-folder"),
  dest = "output_folder",
  action = "store",
  default = getwd(),
  type = "character",
  help = paste0("Folder to save outputs. [default: current working dir]"))
arguments <- parse_args2(parser)
options <- arguments$options
options$positionals <- arguments$args

fun <- switch(options$method,
  avg = Matrix::rowMeans,
  sum = Matrix::rowSums,
  Matrix::rowMeans
)

#TODO: 在 Seurat 应用 MNN 整合数据之间，就查看下两个样本的相关系数。

## Normalization and Variance Stabilization
obj_list <- options$positionals
names(obj_list) <- basename(options$positionals)
obj_list <- lapply(obj_list, function(base_dir) {
  expr_matrix <- read_rhapsody_wta(base_dir, TRUE)
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = expr_matrix, project = basename(base_dir))
  seurat_obj$stim <- basename(base_dir)
  # TODO: 增加 SCTransform 的选项！
  if (isTRUE(options$sctransform)) {
    seurat_obj <- Seurat::SCTransform(
      seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
    return(seurat_obj)
  } else {
    # Perform LogNormalize
    seurat_obj <- Seurat::NormalizeData(seurat_obj)
    seurat_obj <- Seurat::FindVariableFeatures(
      seurat_obj, selection.method = "vst", nfeatures = 2000)
    seurat_obj <- Seurat::ScaleData(seurat_obj, features = rownames(seurat_obj))
    return(seurat_obj)
  }
})

all_stim <- names(obj_list)
comb_list <- as.data.frame(combn(all_stim, 2), stringsAsFactors = FALSE)
names(comb_list) <- NULL
comb_list <- as.list(comb_list)

corr_list <- lapply(comb_list, function(x) {
  stopifnot(length(x) == 2)

  ctrl <- x[[1]]
  test <- x[[2]]

  ctrl_data <- Seurat::GetAssayData(obj_list[[ctrl]])
  test_data <- Seurat::GetAssayData(obj_list[[test]])

  # Gene average expression. (row = genes, col = cells)
  # `expm1` is used in Seurat::AverageExpression
  ctrl_data <- log1p(Matrix::rowMeans(expm1(ctrl_data)))
  test_data <- log1p(Matrix::rowMeans(expm1(test_data)))
  stopifnot(ctrl_data >= 0L)
  stopifnot(test_data >= 0L)

  ctrl_genes <- names(ctrl_data)
  test_genes <- names(test_data)
  all_genes <- union(ctrl_genes, test_genes)

  # Fill NA with zero.
  ctrl_data[setdiff(all_genes, ctrl_genes)] <- 0
  test_data[setdiff(all_genes, test_genes)] <- 0
  stopifnot(length(all_genes) == length(ctrl_data))
  stopifnot(length(all_genes) == length(test_data))

  # Re-order named vector
  ctrl_data <- ctrl_data[all_genes]
  test_data <- test_data[all_genes]
  stopifnot(names(ctrl_data) == names(test_data))

  pearson_r <- cor(ctrl_data, test_data,
    method = "pearson", use = "complete.obs")

  pearson_r^2
})

# Prepare data frame for ggplot2
data_to_plot <- expand.grid(V1 = names(obj_list), V2 = names(obj_list))
coef <- apply(data_to_plot, 1, function(x) {
  names(x) <- NULL
  if (x[[1]] == x[[2]])
    return(1)
  check_result <- lapply(comb_list, function(comb) {
    c(identical(x, comb), identical(rev(x), comb))
  })
  pos <- which(sapply(check_result, any))
  corr_list[[pos]]
})
data_to_plot$coef <- coef

print(data_to_plot)

p <- ggplot2::ggplot(data_to_plot, ggplot2::aes(
      V1, V2, fill = coef)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(
    low = "white", high = "red",
    name = "Pearson R^2") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45, vjust = 1, size = 12, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 12)) +
  ggplot2::coord_fixed() +
  ggplot2::geom_text(ggplot2::aes(
      V1, V2, label = sprintf("%.2f", coef)), color = "black", size = 4) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank())

ggplot2::ggsave(
  file.path(options$output_folder, "corrplot.png"), p,
  width = length(options$positionals),
  height = length(options$positionals))
