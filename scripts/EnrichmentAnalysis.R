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
  "tibble",
  "readr",
  "clusterProfiler",
  "enrichplot",
  "ggplot2",
  "forcats",
  "cowplot",
  "patchwork"
))

library(magrittr)

options(stringsAsFactors = FALSE)
ggplot2::theme_set(cowplot::theme_cowplot())

parse_ratio <- function(ratio) {
  ratio <- sub("^\\s*", "", as.character(ratio))
  ratio <- sub("\\s*$", "", ratio)
  numerator <- as.numeric(sub("/\\d+$", "", ratio))
  denominator <- as.numeric(sub("^\\d+/", "", ratio))
  return(numerator / denominator)
}

go_barplot <- function(resdf, x = "GeneRatio", showCategory = 10, title = "") {
  # Filter empty terms
  resdf <- resdf[!is.na(resdf$Description), ]
  resdf <- resdf[resdf$Count != 0, ]
  # Calculate ratios
  resdf$GeneRatio <- parse_ratio(resdf$GeneRatio)
  resdf$BgRatio <- parse_ratio(resdf$BgRatio)
  # Get top N GO terms
  if (showCategory <= nrow(resdf)) {
      resdf <- resdf[1:showCategory, ]
  }
  # Make top Go terms get higher y values, because of reverse levels
  resdf$Description <- factor(resdf$Description,
                          levels = rev(unique(resdf$Description)))
  # Set plot arguments
  colorBy <- "p.adjust"
  size <- "Count"
  # Set x to user specified col and then 
  resdf <- dplyr::mutate(resdf, x = eval(parse(text = x)))
  # Re-order Go term description by x values, 
  idx <- order(resdf[["x"]], decreasing = TRUE)
  resdf$Description <- factor(
    resdf$Description, levels = rev(unique(resdf$Description[idx])))
  # Plot
  offset <- max(resdf[[x]]) * 0.01
  resdf[["p.adjust.log"]] <- -log10(resdf[["p.adjust"]])
  ggplot2::ggplot(resdf,
      ggplot2::aes_string(x, "Description")) +
      ggplot2::geom_col(
        ggplot2::aes_string(alpha = "p.adjust.log"), fill = "deepskyblue") +
      ggplot2::scale_color_viridis_c() +
      ggplot2::scale_size_continuous(range = c(2, 10)) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_alpha_continuous(name = "-Log10(p.adjust)") +
      ggplot2::ylab(NULL) +
      ggplot2::geom_text(
        ggplot2::aes_string(label = "Description"),
        x = offset, hjust = 0) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank())
}

perform_enrichment_analysis <- function(
  genes, term2gene, term2name, go_data, output_folder) {
  # Perform GO enrichment analysis
  enrich_results <- clusterProfiler::enricher(
    genes, TERM2GENE = term2gene, TERM2NAME = term2name)

  # Add GO_Level to results
  go_level_table <- go_data %>%
    dplyr::group_by(GO_ID) %>%
    dplyr::summarize(Level = unique(GO_Level))
  enrich_results_with_level <- enrich_results %>%
    tibble::as_tibble() %>%
    dplyr::left_join(go_level_table, by = c("ID" = "GO_ID")) %>%
    dplyr::select(ID, Description, Level, dplyr::everything())
  enrich_results_with_level <- enrich_results_with_level %>%
    # RichFactor means DEGs of a GO term divided by All Genes of this GO term.
    dplyr::mutate(
      RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) %>%
    dplyr::select(ID:GeneRatio, BgRatio, RichFactor, dplyr::everything())
    #arrange(desc(geneRatio))

  # Write to tsv
  write.table(
    enrich_results_with_level,
    file = file.path(output_folder, "GO_Enrichment_Results.tsv"),
    sep = "\t",
    row.name = FALSE,
    col.names = TRUE,
    quote = FALSE)

  p_list <- list()
  mf <- enrich_results_with_level %>%
      dplyr::filter(Level == "molecular_function")
  if (nrow(mf) > 0) {
    p_list[[length(p_list) + 1]] <- mf %>% go_barplot(x = "RichFactor",
      showCategory = 20, title = "Molecular Function")
  }

  bp <- enrich_results_with_level %>%
    dplyr::filter(Level == "biological_process")
  if (nrow(bp) > 0) {
    p_list[[length(p_list) + 1]] <- bp %>% go_barplot(x = "RichFactor",
      showCategory = 20, title = "Biological Process")
  }

  cc <- enrich_results_with_level %>%
      dplyr::filter(Level == "cellular_component")
  if (nrow(cc) > 0) {
    p_list[[length(p_list) + 1]] <- cc %>% go_barplot(x = "RichFactor",
      showCategory = 20, title = "Cellular Component")
  }

  p <- patchwork::wrap_plots(p_list)

  ggplot2::ggsave(
    file = file.path(output_folder, "GO_Enrichment_Plot.png"),
    plot = p, width = 7 * length(p_list))
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  options <- list(
    positionals = character(0),
    gene_list_file = NULL,
    go_data_file = NULL,
    ensembldb = NULL,
    group_field = NULL,
    split_by = NULL,
    primary_key = "gene",
    output_folder = getwd()
  )

  optind <- 1
  while (optind <= length(args)) {
    switch(args[optind],
      "--gene-list-file" = {
        # Gene list to conduct enrichment analysis, one gene id per line.
        # Note: you need to discard genes without annotation by yourself.
        # Note: if this argument was provided, positional argument will
        #   be ignored.
        optind <- optind + 1
        options$gene_list_file <- args[optind]
      },
      "--go-data-file" = {
        # GO annotation data.
        optind <- optind + 1
        options$go_data_file <- args[optind]
      },
      "--group-field" = {
        # Which field used to split gene list.
        # Note: work with positional argument.
        optind <- optind + 1
        options$group_field <- args[optind]
      },
      "--split-by" = {
        # Split one table into smaller tables
        optind <- optind + 1
        options$split_by <- args[optind]
      },
      "--primary-key" = {
        # Which field are gene list.
        # Note: work with positional argument
        optind <- optind + 1
        options$primary_key <- args[optind]
      },
      "-C" = {
        # Output folder to place results.
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

  if (is.null(options$go_data_file)) {
    stop("Missing --go-data-file")
  }

  if (is.null(options$gene_list_file) && length(options$positionals) == 0) {
    stop("--gene-list-file or a positional argument is required.")
  }

  # Get GO Data
  # quote="" argument is neccessary for reading the complete table,
  # because it disables quoting.
  go_data <- read.table(
    options$go_data_file, header = TRUE, sep = "\t", quote = "") %>%
    tibble::as_tibble()
  # Discard genes without annotation
  go_data <- go_data %>% dplyr::filter(GO_ID != "")
  # Construct TERM2GENE and TERM2NAME
  # TERM2GENE is a data.frame with first column of term ID and second column of
  # corresponding mapped gene and TERM2NAME is a data.frame with first column of
  # term ID and second column of corresponding term name
  go2gene <- go_data %>% dplyr::select(GO_ID, Gene_ID)
  go2name <- go_data %>% dplyr::select(GO_ID, GO_Name)

  if (!is.null(options$gene_list_file)) {
    genes <- readr::read_lines(options$gene_list_file)
    perform_enrichment_analysis(
      genes, go2gene, go2name, go_data, options$output_folder)
  } else {
    df <- readr::read_csv(
      options$positionals[1], comment = "#")
    if (!is.null(options$group_field)) {
      if (!is.null(options$split_by)) {
        tasks <- split(df, df[[options$split_by]])
      } else  {
        tasks <- list(df)
      }

      lapply(tasks, function(task_df) {
        base::split(task_df, task_df[[options$group_field]]) %>%
          lapply(function(x) {
            output_folder <- file.path(
              ifelse(is.null(options$split_by),
                options$output_folder,
                paste(options$output_folder,
                  unique(x[[options$split_by]]), sep = "_")),
              paste(options$group_field,
                unique(x[[options$group_field]]), sep = "_")
            )
            if (!dir.exists(output_folder))
              dir.create(output_folder, recursive = TRUE)
            tryCatch(
              {
                perform_enrichment_analysis(
                  x[[options$primary_key]],
                  go2gene, go2name, go_data,
                  output_folder
                )
              },
              error = function(e) {
                message(sprintf(
                  "Failed to perform enrichment analysis on %s_%s",
                  options$group_field, unique(x[[options$group_field]])
                ))
              }
            )
          })
      }) %>% invisible()


    } else {
      perform_enrichment_analysis(
        df[[options$primary_key]],
        go2gene, go2name, go_data, options$output_folder)
    }
  }
}