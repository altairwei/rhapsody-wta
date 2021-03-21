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

require_dependencies(c(
  "optparse",
  "readr",
  "dplyr",
  "ggplot2",
  "magrittr",
  "cowplot"
))

library(magrittr)
library(optparse)

ggplot2::theme_set(cowplot::theme_cowplot())

parser <- OptionParser()

parser <- add_option(parser,
  c("-o", "--output-image"),
  dest = "output_image",
  action = "store",
  default = "saturation_curve.png",
  type = "character",
  help = paste0("Output image file. [default: %default]"))
parser <- add_option(parser,
  c("-w", "--which-name"),
  dest = "which_name",
  action = "store",
  default = 0,
  type = "integer",
  help = paste0(
    "Which part of input file path should be assigned to ",
    "the group name of data frame. [default: %default]"))

arguments <- parse_args2(parser)
options <- arguments$options
options$positionals <- arguments$args

if (length(options$positionals) < 1) {
  options$stdin <- TRUE
} else {
  options$stdin <- FALSE
}

filename <- options$output_image

if (length(options$positionals) < 2) {
  # Plot single data frame.
  if (length(options$positionals) < 1) {
    stdin_conn <- file("stdin")
    df <- readr::read_csv(paste(collapse = "\n", readLines(stdin_conn)))
    close(stdin_conn)
  } else {
    input_file <- options$positionals[1]
    df <- readr::read_csv(input_file)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(
      depths, detected_genes, color = origin)) +
    ggplot2::geom_point() +
    ggplot2::geom_line()
} else {
  # Plot multiple data frame in one image.
  input_files <- options$positionals
  df_names <- sapply(input_files, function(x) {
    path_parts <- strsplit(x, .Platform$file.sep)[[1]]
    if (options$which_name == 0 ||
          options$which_name > length(path_parts)) {
      x
    } else {
      path_parts[options$which_name]
    }
  }, USE.NAMES = FALSE)

  names(input_files) <- df_names

  df_list <- lapply(seq_along(input_files), function(index) {
    input_df <- readr::read_csv(input_files[[index]])
    dplyr::mutate(input_df, sample = names(input_files)[[index]])
  })

  df <- dplyr::bind_rows(df_list)

  p <- ggplot2::ggplot(df, ggplot2::aes(
      depths, detected_genes, alpha = origin, color = sample)) +
    ggplot2::geom_line()
}


ggplot2::ggsave(filename, p, device = "png")
