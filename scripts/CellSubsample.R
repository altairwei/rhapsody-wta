#!/usr/bin/env Rscript

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
  "magrittr",
  "Matrix"
))

library(magrittr)
library("optparse")

options(stringsAsFactors = FALSE)

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

parser <- OptionParser()
parser <- add_option(parser, c("-r", "--subsample-ratio"),
  type = "double", help = "Ration to subsample, the value must be less than 1.")
parser <- add_option(parser, c("-C", "--output-folder"),
  default = getwd(), help = "Specify where to save subsampled matrix.")
arguments <- parse_args2(parser)
options <- arguments$options

if (length(arguments$args) < 1) {
  stop("Rhapsody WTA results folder is required.")
}

if (!is.double(options$subsample_ratio)) {
  stop("Subsample ratio must be provided by `-r` or `--subsample-ratio`")
}

if (options$subsample_ratio > 1.0) {
  stop("Subsample ratio must be less than 1.")
}

base_dir <- arguments$args[1]

expr <- read_mtx(base_dir)

message("Subsampling...")
n_cell <- ncol(expr)

# Sampling without putting back
sample_idx <- sample(seq_len(n_cell),
  round(n_cell * options$subsample_ratio), replace = FALSE)

sub_expr <- expr[, sample_idx]

if (!dir.exists(options$output_folder))
  dir.create(options$output_folder)

output_base_name <- sprintf(
  "%s_Subsampled_Expression_Matrix.mtx", basename(base_dir))
mtx_file <- file.path(options$output_folder, output_base_name)
colnames_file <- sprintf("%s.colnames", mtx_file)
rownames_file <- sprintf("%s.rownames", mtx_file)

message(sprintf("Writting %s", mtx_file))
invisible(Matrix::writeMM(sub_expr, file = mtx_file))
write(colnames(sub_expr), colnames_file)
write(rownames(sub_expr), rownames_file)