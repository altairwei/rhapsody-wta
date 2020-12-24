#!/usr/bin/env Rscript

check_deps <- function(x) {
  stopifnot(is.character(x))
  dep_status <- suppressPackageStartupMessages(
    sapply(x, requireNamespace)
  )
  pkg_to_install <- names(dep_status[dep_status == FALSE])
  pkg_to_install
}

require_dependencies <- function(x) {
  pkg_to_install <- check_deps(x)
  if (length(pkg_to_install) > 0) {
    install.packages(pkg_to_install)
  }
}

require_bioconductor <- function(x) {
  pkg_to_install <- check_deps(x)
  if (length(pkg_to_install) > 0) {
    if (!suppressPackageStartupMessages(
      requireNamespace("BiocManager")))
      install.packages("BiocManager")
    BiocManager::install(pkg_to_install)
  }
}

require_dependencies(c(
  "magrittr",
  "dplyr",
  "tibble"
))

require_bioconductor(c(
  "biomaRt",
  "AnnotationDbi"
))

library(magrittr)

options(stringsAsFactors = FALSE)

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  options <- list(
    positionals = character(0),
    host = "plants.ensembl.org",
    mart = "plants_mart",
    dataset = "taestivum_eg_gene",
    output = NULL
  )

  optind <- 1
  while (optind <= length(args)) {
    switch(args[optind],
      "--host" = {
        optind <- optind + 1
        options$host <- args[optind]
      },
      "--mart" = {
        optind <- optind + 1
        options$host <- args[optind]
      },
      "--dataset" = {
        optind <- optind + 1
        options$host <- args[optind]
      },
      "-o" = {
        optind <- optind + 1
        options$output <- args[optind]
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

  # Ger all params
  dest_attrs <- c(
    "ensembl_gene_id", "go_id", "name_1006", "namespace_1003"
  )

  # Build connections
  dataset_conn <- biomaRt::useMart(
    options$mart, dataset = options$dataset, host = options$host)

  # Get all annotated gene list
  all_gene_list <- biomaRt::getBM(
    attributes = c("ensembl_gene_id"),
    filters = "chromosome_name",
    values = AnnotationDbi::keys(dataset_conn, keytype = "chromosome_name"),
    mart = dataset_conn
  ) %>% tibble::as_tibble()

  all_gene_go_data <- biomaRt::getBM(
    attributes = dest_attrs,
    filters = "ensembl_gene_id",
    values = all_gene_list[["ensembl_gene_id"]],
    mart = dataset_conn
  ) %>% tibble::as_tibble()

  cat(paste(
    "Total genes: ",
    length(unique(all_gene_go_data[["ensembl_gene_id"]])), "\n"),
    file = stderr())

  all_gene_go_data <- all_gene_go_data %>% dplyr::rename(
    Gene_ID = ensembl_gene_id,
    GO_ID = go_id,
    GO_Name = name_1006,
    GO_Level = namespace_1003
  )

  output <- if (is.null(options$output)) stdout() else options$output

  write.table(
    all_gene_go_data, file = output,
    sep = "\t", row.name = FALSE, col.names = TRUE, quote = FALSE)
}