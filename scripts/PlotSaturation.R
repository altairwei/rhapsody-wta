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
  "readr",
  "ggplot2",
  "magrittr",
  "cowplot"
))

library(magrittr)

ggplot2::theme_set(cowplot::theme_cowplot())

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("ERROR: Image filename is required.")
}

filename <- args[1]

stdin_conn <- file("stdin")
df <- readr::read_csv(paste(collapse = "\n", readLines(stdin_conn)))
close(stdin_conn)

p <- ggplot2::ggplot(df, ggplot2::aes(
    depths, detected_genes, color = origin)) +
  ggplot2::geom_point() +
  ggplot2::geom_line()


ggplot2::ggsave(filename, p, device = "png")
