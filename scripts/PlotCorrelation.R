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
  "magrittr",
  "Matrix",
  "pheatmap",
  "reshape2",
  "rhapsodykit",
  "ComplexHeatmap"
))

library(magrittr)
library(optparse)

options(stringsAsFactors = FALSE)

parser <- OptionParser()
parser <- add_option(parser,
  c("--sum"),
  dest = "average",
  action = "store_false",
  help = "Create pseudo-bulk data by summing the counts from cells.")
parser <- add_option(parser,
  c("--avg"),
  dest = "average",
  action = "store_true",
  default = TRUE,
  help = "Create pseudo-bulk data by averaging the counts from cells.")

parser <- add_option(parser,
  c("--cor-method"),
  dest = "cor_method",
  action = "store",
  default = "pearson",
  help = paste0("Which method should be used to compute ",
    "correlation coefficient (or covariance). [default: %default]"))
parser <- add_option(parser,
  c("--pearson"),
  dest = "cor_method",
  action = "callback",
  callback = function(...) "pearson",
  help = "Compute Pearson's correlation coefficient (or covariance)")
parser <- add_option(parser,
  c("--kendall"),
  dest = "cor_method",
  action = "callback",
  callback = function(...) "kendall",
  help = "Compute Kendall's correlation coefficient (or covariance)")
parser <- add_option(parser,
  c("--spearman"),
  dest = "cor_method",
  action = "callback",
  callback = function(...) "spearman",
  help = "Compute Spearman's correlation coefficient (or covariance)")

parser <- add_option(parser,
  c("--normalization-method"),
  dest = "norm_method",
  action = "store",
  default = "LogNormalize",
  help = paste0("Specify the method used to normalize expression ",
    "data. [default: %default]"))
parser <- add_option(parser,
  c("--log"),
  dest = "norm_method",
  action = "callback",
  callback = function(...) "LogNormalize",
  help = paste0("Use LogNormalize of Seurat::NormalizeData ",
    "to perform normalization."))
parser <- add_option(parser,
  c("--clr"),
  dest = "norm_method",
  action = "callback",
  callback = function(...) "CLR",
  help = paste0("Use CLR of Seurat::NormalizeData ",
    "to perform normalization."))
parser <- add_option(parser,
  c("--rc"),
  dest = "norm_method",
  action = "callback",
  callback = function(...) "RC",
  help = paste0("Use RC of Seurat::NormalizeData ",
    "to perform normalization."))
parser <- add_option(parser,
  c("--sct"),
  dest = "norm_method",
  action = "callback",
  callback = function(...) "SCTransform",
  help = "Use Seurat::SCTransform to perform normalization.")

parser <- add_option(parser,
  c("-O", "--output-folder"),
  dest = "output_folder",
  action = "store",
  default = getwd(),
  type = "character",
  help = paste0("Folder to save outputs. [default: current working dir]"))
parser <- add_option(parser,
  c("-o", "--name-prefix"),
  dest = "name_prefix",
  action = "store",
  default = "Sample_Correlation",
  type = "character",
  help = "The name prefix of output file. [default: %default]")
arguments <- parse_args2(parser)
options <- arguments$options
options$positionals <- arguments$args

fun_to_apply <- if (options$average) Matrix::rowMeans else Matrix::rowSums

#TODO: 对相关性热图进行聚类分析。
#TODO: 增加选项输出 psuedo-bulk RNA-Seq 数据

expr_df <- rhapsodykit::make_psuedo_bulk(
  options$positionals,
  normalization = options$norm_method,
  method = ifelse(options$average, "avg", "sum"))

cormat <- cor(expr_df, method = options$cor_method, use = "complete.obs")

data_to_plot <- reshape2::melt(cormat, value.name = "coef")

write.table(
  data_to_plot, file = file.path(
    options$output_folder,
    sprintf("%s.tsv", options$name_prefix)),
  sep = "\t", row.name = FALSE, col.names = TRUE, quote = FALSE)

scale_fill_name <- switch(options$cor_method,
  pearson = "Pearson's r",
  kendall = "Kendall's tau",
  spearman = "Spearman's rho"
)

# p <- ggplot2::ggplot(data_to_plot, ggplot2::aes(
#       Var1, Var2, fill = coef)) +
#   ggplot2::geom_tile() +
#   ggplot2::scale_fill_gradient(
#     low = "white", high = "red",
#     name = scale_fill_name) +
#   ggplot2::theme_minimal() +
#   ggplot2::theme(
#     axis.text.x = ggplot2::element_text(
#       angle = 45, vjust = 1, size = 12, hjust = 1),
#     axis.text.y = ggplot2::element_text(size = 12)) +
#   ggplot2::coord_fixed() +
#   ggplot2::geom_text(ggplot2::aes(
#       Var1, Var2, label = sprintf("%.2f", coef)), color = "black", size = 4) +
#   ggplot2::theme(
#     axis.title.x = ggplot2::element_blank(),
#     axis.title.y = ggplot2::element_blank(),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.border = ggplot2::element_blank(),
#     panel.background = ggplot2::element_blank(),
#     axis.ticks = ggplot2::element_blank())

# ggplot2::ggsave(
#   file.path(
#     options$output_folder,
#     sprintf("%s.png", options$name_prefix)),
#   p,
#   width = length(options$positionals),
#   height = length(options$positionals))

p <- ComplexHeatmap::Heatmap(
  cormat,
  name = scale_fill_name,
  col = c("white", "red"),
  row_names_side = "left",
  row_dend_side = "left",
  column_names_side = "bottom",
  column_dend_side = "top",
  column_names_rot = 45,
  show_column_dend = FALSE,

  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(
      sprintf("%.2f", cormat[i, j]),
      x, y, gp = grid::gpar(fontsize = 14)
    )
  }
)

pngfile <- file.path(
  options$output_folder, sprintf("%s.png", options$name_prefix))
png(filename = pngfile,
  width = length(options$positionals),
  height = length(options$positionals),
  units = "in",
  res = 300)
ComplexHeatmap::draw(p)
invisible(dev.off())