#!/usr/bin/env Rscript

library(optparse)
library(Seurat)

parser <- OptionParser()

parser <- add_option(parser,
  c("--group-rule"),
  dest = "group_rule",
  action = "store",
  default = NULL,
  type = "character",
  help = paste0("The pattern to extract group information from sample name.")
)

arguments <- parse_args2(parser)
options <- arguments$options
options$positionals <- arguments$args

if (length(options$positionals) < 1) {
  stop("At least one position argument is required.\n")
}

for (rds_filename in options$positionals) {
  obj <- readRDS(rds_filename)

  stopifnot(inherits(obj, "Seurat"))

  # Replace outdated keys
  if (!is.null(obj@meta.data$stim)) {
    obj$sample <- obj$stim
    obj$stim <- NULL
  }

  # Extract group
  if (!is.null(options$group_rule)) {
    obj$group <- stringr::str_extract(
      obj$sample, options$group_rule)
  }

  saveRDS(obj, paste0(rds_filename, ".new"))

  gc()
}