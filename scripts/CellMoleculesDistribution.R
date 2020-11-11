#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

options <- list(
  flags = args[grep("^-", args)],
  positionals = args[grep("^-", args, invert = TRUE)]
)

if (length(options$positionals) != 2) {
  stop("Two position arguments are required.\n")
}

options$data_file <- options$positionals[1]
options$output_file <- options$positionals[2]

suppressMessages(library(tidyverse))

df <- readr::read_tsv(options$data_file, comment = "#",
  col_types = cols(
    Cell_Index = col_factor(),
    Gene = col_character(),
    RSEC_Reads = col_double(),
    Raw_Molecules = col_double(),
    RSEC_Adjusted_Molecules = col_double()
  )
)

cell_df <- df %>%
  group_by(Cell_Index) %>%
  summarise(Molecules = sum(RSEC_Adjusted_Molecules))

p <- cell_df %>%
  ggplot(aes(Molecules)) +
  geom_histogram(bins = 100)

ggsave(options$output_file, p)