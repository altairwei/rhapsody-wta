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
  "magrittr",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "readr",
  "patchwork"
))

library(magrittr)

options(stringsAsFactors = FALSE)
ggplot2::theme_set(cowplot::theme_cowplot())

save_plot <- function(filename, plot, width = 7, height = 7, ...) {
  message(sprintf("Saving %d x %d in image: %s", width, height, filename))
  ggplot2::ggsave(filename, plot, width = width, height = height, ...)
}

search_file <- function(file_pattern, base_dir) {
  if (!dir.exists(base_dir)) {
    stop("Directory provided does not exist")
  }

  filename <- Sys.glob(
    file.path(base_dir, file_pattern))

  if (length(filename) != 1) {
    stop(sprintf("`%s` missing or more than one file.", file_pattern))
  }

  filename
}

read_annot_mol <- function(base_dir = ".") {
  annot_mol_file <- search_file("*_Annotation_Molecule.csv.gz", base_dir)
  message("Reading ", annot_mol_file)
  readr::read_csv(
    gzfile(annot_mol_file), progress = TRUE,
    col_names = c(
      "cell", "umi", "gene", "umi_count", "rsec_corrected_umi_count",
      "dbec_corrected_umi_count", "parent_umi"),
    col_types = readr::cols(
      cell = readr::col_character(),
      umi = readr::col_character(),
      gene = readr::col_character(),
      umi_count = readr::col_double(),
      rsec_corrected_umi_count = readr::col_double(),
      dbec_corrected_umi_count = readr::col_double(),
      parent_umi = readr::col_character()
    )
  )
}

read_cell_index <- function(base_dir = ".") {
  putative_cells_file <- search_file("*_Putative_Cells_Origin.csv", base_dir)
  message("Reading ", putative_cells_file)
  readr::read_csv(
    putative_cells_file, comment = "#",
    col_types = readr::cols(
      Cell_Index = readr::col_character(),
      Algorithm = readr::col_character()
    )
  )
}

read_umi_stats <- function(base_dir = ".") {
  umi_stats_file <- search_file("*_UMI_Adjusted_Stats.csv", base_dir)
  message("Reading ", umi_stats_file)
  readr::read_csv(
    umi_stats_file, comment = "#",
    col_types = readr::cols(
      Gene = readr::col_character(),
      Status = readr::col_character(),
      Raw_Reads = readr::col_integer(),
      Raw_Molecules = readr::col_integer(),
      Raw_Seq_Depth = readr::col_double(),
      RSEC_Adjusted_Molecules = readr::col_integer(),
      RSEC_Adjusted_Seq_Depth = readr::col_double(),
      RSEC_Adjusted_Seq_Depth_without_Singletons = readr::col_double()
    )
  )
}

read_expr_data <- function(base_dir = ".") {
  expressiong_data_file <- search_file("*_Expression_Data.st", base_dir)
  message("Reading ", expressiong_data_file)
  readr::read_tsv(expressiong_data_file, comment = "#",
    col_types = readr::cols(
      Cell_Index = readr::col_character(),
      Gene = readr::col_character(),
      RSEC_Reads = readr::col_double(),
      Raw_Molecules = readr::col_double(),
      RSEC_Adjusted_Molecules = readr::col_double()
    )
  )
}

umi_histogram <- function(data) {
  raw_umi <- nrow(data)
  rsec_umi <- sum(data$rsec_corrected_umi_count != 0)
  line1 <- sprintf("Raw UMI Count: %d", raw_umi)
  line2 <- sprintf("RSEC UMI Count: %d (%.2f%%)",
    rsec_umi, (rsec_umi / raw_umi) * 100)
  df <- dplyr::filter(data, rsec_corrected_umi_count != 0)
  total_num <- nrow(df)
  depth_gt_1 <- nrow(
    df <- dplyr::filter(df, rsec_corrected_umi_count > 1)
  )
  line3 <- sprintf(
    "Depth > 1: %.2f%%",
    round((depth_gt_1 / total_num) * 100, digits = 2))
  percent_label <- paste(line1, line2, line3, sep = "\n")
  ggplot2::ggplot(df, ggplot2::aes(rsec_corrected_umi_count)) +
    ggplot2::geom_histogram(binwidth = 1) +
    ggplot2::scale_x_continuous(limits = c(1, 40)) +
    ggplot2::annotate(
      geom = "text", label = percent_label, size = 4,
      x = 20, y = Inf, hjust = "left", vjust = "top") +
    ggplot2::xlab("Sequencing Depth of RSEC Molecules") +
    ggplot2::ylab("Number of Molecules")
}

cell_umi_saturation_histogram <- function(mol_df, cell_df) {
  df <- dplyr::filter(mol_df, rsec_corrected_umi_count != 0)
  compute_pct <- function(x) {
    stopifnot(is.numeric(x))
    total_num <- length(x)
    depth_gt_1 <- sum(x > 1)
    percent <- (depth_gt_1 / total_num) * 100
    percent
  }
  percent_df <- df %>%
    dplyr::group_by(cell) %>%
    dplyr::summarise(pct = compute_pct(rsec_corrected_umi_count))
  percent_df <- percent_df[percent_df$cell %in% cell_df$Cell_Index, ]
  percent_df <- dplyr::left_join(
    percent_df, cell_df, by = c("cell" = "Cell_Index"))

  ggplot2::ggplot(percent_df, ggplot2::aes(pct, fill = Algorithm)) +
    ggplot2::geom_histogram(position = "stack", bins = 100) +
    ggplot2::scale_x_continuous(limits = c(0, 100)) +
    ggplot2::ggtitle("Cell UMIs Saturation Distribution") +
    ggplot2::xlab("Percentage of RSEC molecules which depth > 1") +
    ggplot2::ylab("Count of Cells") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

gene_seq_depth_histogram <- function(umi_stats_df) {
  p_seq_depth <- ggplot2::ggplot(umi_stats_df,
    ggplot2::aes(RSEC_Adjusted_Seq_Depth, fill = Status)) +
    ggplot2::geom_histogram(bins = 80) +
    ggplot2::scale_x_continuous(limits = c(0, 30)) +
    ggplot2::ggtitle("Gene Sequence Depth Distribution") +
    ggplot2::xlab("RSEC Adjusted Seq Depth") +
    ggplot2::ylab("Number of Genes") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  p_umi_count <- ggplot2::ggplot(umi_stats_df,
    ggplot2::aes(RSEC_Adjusted_Molecules, fill = Status)) +
    ggplot2::geom_histogram(bins = 80) +
    ggplot2::scale_x_continuous(limits = c(0, 2500)) +
    ggplot2::ggtitle("Gene UMIs Count Distribution") +
    ggplot2::xlab("Count of RSEC Adjusted Molecules") +
    ggplot2::ylab("Number of Genes") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  patchwork::wrap_plots(list(p_seq_depth, p_umi_count)) +
    patchwork::plot_layout(guides = "collect")
}

cell_gene_count_histogram <- function(expr_df) {
  cell_df <- expr_df %>%
    dplyr::group_by(Cell_Index) %>%
    dplyr::summarise(
      Molecules = sum(RSEC_Adjusted_Molecules),
      GeneCount = length(unique(Gene)))

  p_cell_gene_count <- ggplot2::ggplot(cell_df,
    ggplot2::aes(GeneCount)) +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::scale_x_continuous(limits = c(0, 10000)) +
    ggplot2::ggtitle("Cell Detected Gene Number Distribution") +
    ggplot2::xlab("Detected Gene Number") +
    ggplot2::ylab("Number of Cells") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  p_cell_umi_count <- ggplot2::ggplot(cell_df,
    ggplot2::aes(Molecules)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::scale_x_continuous(limits = c(0, 40000)) +
    ggplot2::ggtitle("Cell UMIs Count Distribution") +
    ggplot2::xlab("Number of Molecules") +
    ggplot2::ylab("Number of Cells") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  patchwork::wrap_plots(list(p_cell_gene_count, p_cell_umi_count))
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  options <- list(
    positionals = character(0),
    process = 1,
    output_folder = getwd()
  )

  optind <- 1
  while (optind <= length(args)) {
    switch(args[optind],
      "--process" = {
        optind <- optind + 1
        options$process <- as.integer(args[optind])
      },
      "-C" = {
        optind <- optind + 1
        options$output_folder <- args[optind]
      },
      {
        if (startsWith(args[optind], "-")) {
          stop(sprintf("Unknown option: %s", args[optind]))
        } else {
          # CSV table contain gene list and group info
          # Note: you should filter out un-annotated genes by yourself.
          options$positionals <- append(options$positionals, args[optind])
        }
      }
    )
    optind <- optind + 1
  }

  if (length(options$positionals) < 1) {
    stop("At least one position argument is required.\n")
  }

  if (options$process > 1) {
    if (suppressPackageStartupMessages(!requireNamespace("future")))
      install.packages("future")
    future::plan("multiprocess", workers = options$process)
    # Set global size to 2GB
    options(future.globals.maxSize = 4 * 1024^3)
  }

  invisible(
    lapply(options$positionals, function(base_folder) {
      if (!dir.exists(options$output_folder))
        dir.create(options$output_folder, recursive = TRUE)
      tryCatch(
        {
          annot_mol_df <- read_annot_mol(base_folder)
          cell_index_df <- read_cell_index(base_folder)
          p_umi_hist <- umi_histogram(annot_mol_df)
          save_plot(
            file.path(options$output_folder, "QC_Molecules_Depth.png"),
            p_umi_hist
          )
          p_umi_sat <- cell_umi_saturation_histogram(
            annot_mol_df, cell_index_df)
          save_plot(
            file.path(options$output_folder, "QC_Cell_UMI_Saturation.png"),
            p_umi_sat
          )
        },
        error = function(e) {
          message(toString(e))
        }
      )

      tryCatch(
        {
          umi_stats_df <- read_umi_stats(base_folder)
          p_gene_seq_depth <- gene_seq_depth_histogram(umi_stats_df)
          save_plot(
            file.path(options$output_folder, "QC_Gene_Seq_Depth.png"),
            p_gene_seq_depth, width = 14
          )
        },
        error = function(e) {
          message(toString(e))
        }
      )

      tryCatch(
        {
          expr_data_df <- read_expr_data(base_folder)
          p_cell_gene_count <- cell_gene_count_histogram(expr_data_df)
          save_plot(
            file.path(options$output_folder, "QC_Cell_Gene_Distribution.png"),
            p_cell_gene_count, width = 14
          )
        },
        error = function(e) {
          message(toString(e))
        }
      )
    })
  )

}