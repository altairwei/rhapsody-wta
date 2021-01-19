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
  "magrittr",
  "Matrix",
  "Seurat",
  "sctransform",
  "devtools"
))

library(magrittr)
library(optparse)

load_sources()

options(stringsAsFactors = FALSE)

parser <- OptionParser()
parser <- add_option(parser, c("-m", "--method"), dest = "method",
  action = "store", default = "avg", type = "character",
  help = paste0("Create pseudo-bulk data by summing/averaging",
    " the counts from cells. Choices are `avg` and `sum`. [default: %default]"))
parser <- add_option(parser, c("--sum"), dest = "method",
  action = "callback", callback = function(...) "sum",
  help = "Create pseudo-bulk data by summing the counts from cells.")
parser <- add_option(parser, c("--avg"), dest = "method",
  action = "callback", callback = function(...) "avg",
  help = "Create pseudo-bulk data by averaging the counts from cells.")
parser <- add_option(parser, c("--sctransform"), action = "store_true",
  help = "Use Seurat::SCTransform to perform normalization.")
arguments <- parse_args2(parser)
options <- arguments$options
positionals <- arguments$args

fun <- switch(options$method,
  avg = rowMeans, sum = rowSums, rowMeans
)

#TODO: 在 Seurat 应用 MNN 整合数据之间，就查看下两个样本的相关系数。

## Normalization and Variance Stabilization
print(read_mtx)