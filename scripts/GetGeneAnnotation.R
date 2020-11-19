#!/usr/bin/env Rscript

suppressMessages(library(biomaRt))
suppressMessages(library(tibble))
suppressMessages(library(readr))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))

table_collapse <- function(x) {
    empty <- ""
    if (all(x == empty)) {
        return("-")
    } else {
        return(paste(unique(x), collapse = "|"))
    }
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  options <- list(
    positionals = character(0),
    host = "plants.ensembl.org",
    mart = "plants_mart",
    dataset = "taestivum_eg_gene",
    primary_key = "gene",
    inplace = FALSE
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
      "--primary-key" = {
        optind <- optind + 1
        options$host <- args[optind]
      },
      "--inplace" = {
        options$inplace <- TRUE
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

  df <- readr::read_csv(
    options$positionals[1], comment = "#", progress = TRUE)

  dataset_conn <- biomaRt::useMart(
    options$mart, dataset = options$dataset, host = options$host)
  annot_df <- biomaRt::getBM(
      attributes = c(
        "ensembl_gene_id", "description", "go_id", "name_1006",
        "namespace_1003", "plant_reactome_pathway", "plant_reactome_reaction",
        "interpro_short_description"),
      filters = "ensembl_gene_id",
      values = unique(df[[options$primary_key]]),
      mart = dataset_conn
  )

  annot_df <- annot_df %>%
    tibble::as_tibble() %>%
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::summarize(
      Description = table_collapse(description),
      GO_ID = table_collapse(go_id),
      GO_Term = table_collapse(name_1006),
      GO_Domain = table_collapse(namespace_1003),
      Interpro_Short_Description = table_collapse(interpro_short_description),
      Reactome_Pathway = table_collapse(plant_reactome_pathway),
      Reactome_Reaction = table_collapse(plant_reactome_reaction)
  )

  output <- dplyr::left_join(
    df, annot_df, by = c("gene" = "ensembl_gene_id")) %>%
    readr::format_csv()

  writeLines(output, ifelse(isTRUE(options$inplace),
    options$positionals[1], stdout()))
}