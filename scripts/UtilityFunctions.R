# Make {renv} records pkg {R.utils}
loadNamespace("R.utils")

#' Display cluster markers in a html table
#'
#' @param df Data frame returned by Seurat::FindAllMarkers
marker_table <- function(df) {
  df %>%
    dplyr::select(gene, dplyr::everything()) %>%
    reactable::reactable(
      defaultColDef = reactable::colDef(
        minWidth = 40
      ),
      columns = list(
        gene = reactable::colDef(minWidth = 80),
        p_val = reactable::colDef(
          cell = function(x) format(x, digits=3, scientific = TRUE)
        ),
        p_val_adj = reactable::colDef(
          cell = function(x) format(x, digits=3, scientific = TRUE)
        ),
        avg_log2FC = reactable::colDef(
          format = reactable::colFormat(
            digits = 2
          )
        ),
        pct.1 = reactable::colDef(
          format = reactable::colFormat(
            percent = TRUE,
            digits = 1
          )
        ),
        pct.2 = reactable::colDef(
          format = reactable::colFormat(
            percent = TRUE,
            digits = 1
          )
        ),
        .rownames = reactable::colDef(
          show = FALSE
        )
      )
    )
}

plot_bootstrap_distribution <- function(
  res, clusters = NULL,
  facet_by = "cellTypes", ncol = NULL, nrow = NULL
) {
  # Use BCa as CI
  df_to_plot <- res$results %>%
    tibble::as_tibble() %>%
    dplyr::filter(method == "BCa") %>%
    dplyr::mutate(
      time = sapply(strsplit(cond, split = "-"), "[", 1),
      rep = sapply(strsplit(as.character(subject), split = "-"),
                   function(x) paste(x[2:3], collapse = "-")),
      treatment = sapply(strsplit(cond, split = "-"), "[", 2)
    )
  
  prop_df <- res$thetastar %>%
    as.data.frame()
  
  colnames(prop_df) <- paste("BS_", seq_len(ncol(res$thetastar)), sep = "")
  
  df_to_plot <- dplyr::bind_cols(df_to_plot, prop_df) %>%
    tidyr::pivot_longer(tidyr::starts_with("BS_"), names_to = "bootstrap", values_to = "prop")
  
  #df_to_plot <- df_to_plot %>%
  #  dplyr::group_by(time, treatment, cellTypes) %>%
  #  dplyr::summarise(mean = mean(prop), sd = sd(prop))
  
  df_to_plot <- dplyr::bind_rows(
    df_to_plot,
    df_to_plot %>%
      dplyr::filter(time == "0DPI") %>%
      dplyr::mutate(treatment = "PNR2"),
    df_to_plot %>%
      dplyr::filter(time == "0DPI") %>%
      dplyr::mutate(treatment = "TR4")
  )
  
  if (!is.null(clusters)) {
    df_to_plot <- dplyr::filter(df_to_plot, cellTypes %in% clusters)
    df_to_plot$cellTypes <- factor(df_to_plot$cellTypes, levels = clusters)
  }
  
  p <- ggplot2::ggplot(df_to_plot,
                       ggplot2::aes(x = time, y = prop, color = treatment)) +
    ggplot2::geom_boxplot(
      width = 0.2, outlier.shape = NA,
      position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(group = treatment),
      stat = "summary", fun = median,
      position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_bw() +
    NULL
  
  if (!is.null(facet_by))
    p <- p + ggplot2::facet_wrap(facet_by, scales = "free_y", ncol = ncol, nrow = nrow, drop = FALSE)
  
  p
}

plot_markers <- function(obj, ...) {
  Seurat::DotPlot(
    obj,
    features = list(
      # Epidermal Cells
      "FDH" = c(
        "TraesCS4B02G297500",
        "TraesCS4D02G296400",
        "TraesCS4A02G007400"
      ),
      "ATML1" = c(
        "TraesCS2A02G474000",
        "TraesCS2D02G473700"
      ),
      "DCR" = c(
        "TraesCS1A02G341300",
        "TraesCS1D02G343400"
      ),
      
      # EP3 是排水孔相关基因
      "EP3" = c(
        "TraesCS2A02G350700",
        "TraesCS2D02G348800",
        "TraesCS6D02G199500",
        "TraesCS6A02G216100"
      ),
      
      # Guardian Cells
      "ALMT12" = c(
        "TraesCS1D02G194000",
        "TraesCS1A02G189900",
        "TraesCS1B02G192000"
      ),
      "MYB60" = c(
        "TraesCS4A02G322200",
        "TraesCS5D02G552200"
      ),
      "HIC" = c(
        "TraesCS4D02G226100"
      ),
      
      # Mesophyll Cells
      "RBCS" = c(
        "TraesCS2A02G066800",
        #"TraesCS2B02G079200",
        #"TraesCS2D02G065200",
        #"TraesCS2D02G065300",
        "TraesCS5A02G165400",
        #"TraesCS5A02G165700",
        #"TraesCS5B02G162600",
        #"TraesCS5B02G162800",
        #"TraesCS5D02G169600",
        "TraesCS5D02G169900"
      ),
      "CAB3" = c(
        "TraesCS7A02G276400",
        #"TraesCS1D02G411300",
        "TraesCS1B02G317500",
        #"TraesCS7D02G276300",
        #"TraesCS5B02G353200",
        #"TraesCS5A02G350600",
        "TraesCS1A02G403300"
      ),
      "LHCB2.1" = c(
        "TraesCS5D02G329200",
        "TraesCS5B02G322900",
        "TraesCS5A02G322500"
      ),
      "CA1" = c(
        #"TraesCS7D02G443400",
        "TraesCS7B02G354800",
        "TraesCS3A02G230000",
        #"TraesCS3D02G223300",
        "TraesCS3B02G259300"
      ),
      "AOC2" = c(
        "TraesCS6D02G314300",
        "TraesCS6A02G334800",
        "TraesCS6B02G365200"
      ),
      
      # Vascular Cells
      "SULTR3;4" = c(
        "TraesCS7A02G088700",
        "TraesCS4A02G388000",
        "TraesCS7D02G084100"
      ),
      "TaGSr" = c(
        "TraesCS4B02G240900",
        "TraesCS4D02G240700",
        "TraesCS4A02G063800"
      ),
      "gl-OXO" = c(
        "TraesCS4D02G032000",
        "TraesCS4B02G033300",
        "TraesCS4A02G279200",
        "TraesCS4D02G031800"
      ),
      "TaSUT1" = c(
        "TraesCS4A02G016400",
        "TraesCS4B02G287800",
        "TraesCS4D02G286500"
      ),
      "CPIII" = c(
        "TraesCS6B02G050700",
        "TraesCS6D02G041700",
        "TraesCS6A02G036100"
      )
    ),
    ...
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text.x = ggplot2::element_text(angle = 90)
    )
}

choose_indent_res <- function(obj, meta_col) {
  Seurat::Idents(obj) <- meta_col
  lvs <- type.convert(levels(obj))
  Seurat::Idents(obj) <- factor(Seurat::Idents(obj), levels = sort(lvs))
  obj
}

centerPopulationStimData <- function(sce) {
  sce <- sce[, sce$cluster_id %in% c(
    paste0("Me_", 1:6), paste0("Va_", 1:4),
    "MPV_1", "MPV_2", "BS")]

  embed_umap <- reducedDim(sce, "UMAP")

  sce <- sce[,
     embed_umap[, 1L] < 8
     & embed_umap[, 1L] > -8
     & embed_umap[, 2L] < 8
  ]

  sce$cluster_id <- droplevels(sce$cluster_id)
  if (!is.null(sce$cellType))
    sce$cellType <- droplevels(sce$cellType)

  sce
}

centerPopulationMockData <- function(sce) {
  sce <- sce[, sce$cluster_id %in% c(
    "Me_α", "Me_β", "Me_γ", "Me_δ", "Me_ε",
    "BS", "MPV", "Va_α", "Va_β", "Va_γ", "Va_δ"
  )]

  embed_umap <- reducedDim(sce, "UMAP")

  sce <- sce[,
    embed_umap[, 1L] < 10
    & embed_umap[, 1L] > -10
    & embed_umap[, 2L] < 10
  ]
  
  sce$cluster_id <- droplevels(sce$cluster_id)
  if (!is.null(sce$cellType))
    sce$cellType <- droplevels(sce$cellType)

  sce
}

#' Deconvolution of LCM Samples
#'
#' @param seurat Seurat object
#' @param lcm_file Count file of LCM samples
#' @param markers Marker genes used for deconvolution
#' @param design Named vector to specify sample information, such
#' as which duplicates belong to the same group (case/control)
#' @param ... Arguments passed to \link[RNAMagnet]{runCIBERSORT}
#'
#' @return A data frame returned by \link[RNAMagnet]{runCIBERSORT}
deconvLCM <- function(seurat, lcm_file, markers, design, ...) {
  # Load LCM data
  lcm_counts <- readr::read_tsv(lcm_file, comment = "#")
  lcm_counts <- lcm_counts[!(names(lcm_counts) %in% c("Chr", "Start", "End", "Strand", "Length"))]
  names(lcm_counts) <- c("Gene", basename(dirname(names(lcm_counts)[-1])))
  lcm_mtx <- as.matrix(lcm_counts[-1])
  rownames(lcm_mtx) <- lcm_counts$Gene
  lcm_mtx <- lcm_mtx[ rowSums(lcm_mtx) > 1, ]

  # scRNA-seq data
  # 是否需要使用归一化后的表达之呢？不需要，MuSiC 就是 raw counts
  sc_counts <- SeuratObject::GetAssayData(seurat, assay = "RNA", slot = "counts")
  sc_mtx <- vapply(
    X = levels(Seurat::Idents(seurat)),
    FUN = function(i) {
      Matrix::rowMeans(
        sc_counts[markers, Seurat::WhichCells(seurat, ident = i)])
    },
    FUN.VALUE = numeric(length = length(markers))
  )
  
  sc_mtx <- sc_mtx[ rowSums(sc_mtx) > 1, ]

  CIBER <- RNAMagnet::runCIBERSORT(
    exprs = lcm_mtx,
    base = sc_mtx,
    design = LCM_design,
    markergenes = intersect(rownames(lcm_mtx), rownames(sc_mtx)),
    ...
  )

  CIBER
}

#' Scatter Plot of Cell Type Estimate
#'
#' @inheritParams deconvLCM
#'
#' @return ggplot object
deconvScatter <- function(CIBER, design, nrow = NULL, ncol = NULL, theme_size = 12) {
  ggplot2::ggplot(
    data = CIBER,
    mapping = ggplot2::aes(
      x = factor(SampleClass, unique(design)),
      y= Fraction, color = as.character(CellType))
  ) +
    ggplot2::geom_point(stat = "summary", fun = mean) +
    ggplot2::geom_errorbar(
      stat="summary",
      fun.min = function(x) mean(x)+sd(x)/sqrt(length(x)),
      fun.max = function(x) mean(x)-sd(x)/sqrt(length(x)),
      width = 0.2
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::facet_wrap(~ CellType, scales = "free_y",
                        nrow = nrow, ncol = ncol) +
    ggplot2::theme_bw(base_size=theme_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::guides(color = "none") +
    ggplot2::ylab("CIBERSORT estimate (a.u.)") +
    ggplot2::xlab("Niche")
}

#' Heatmap of Cell Type Estimate
#'
#' @inheritParams deconvLCM
#'
#' @return A \link[ComplexHeatmap]{Heatmap} object
deconvHeatmap <- function(CIBER) {
  mtx_df <- CIBER %>%
    dplyr::group_by(SampleClass, CellType) %>%
    dplyr::summarise(
      MeanFrac = mean(Fraction)
    ) %>%
    dplyr::select(CellType, SampleClass, MeanFrac) %>%
    tidyr::pivot_wider(names_from = SampleClass, values_from = MeanFrac)

  mtx <- as.matrix(mtx_df[, -1])
  rownames(mtx) <- mtx_df$CellType

  ComplexHeatmap::Heatmap(
    matrix = mtx,
    col = c("white", "blue", "red"),
    row_dend_width = grid::unit(4, "cm"),
    column_order = unique(LCM_design),
    heatmap_legend_param = list(
      title = "Fraction",
      #title_position = "leftcenter",
      legend_direction = "vertical",
      legend_height = grid::unit(4, "cm")
    )
  )
}


pseudobulk <- function(seurat,
    type = c("counts", "logcounts", "cpm", "vstresiduals"),
    fun = NULL) {
  seurat$cellType <- Seurat::Idents(seurat)
  sce <- Seurat::as.SingleCellExperiment(seurat, assay = "RNA")
  sce <- muscat::prepSCE(
    sce,
    kid = "ident", # subpopulation assignments
    gid = "group", # group IDs (ctrl/stim)
    sid = "sample", # sample IDs (ctrl/stim.1234)
    drop = FALSE
  )

  sce <- sce[Matrix::rowSums(SingleCellExperiment::counts(sce) > 0) > 0, ]

  type <- match.arg(type)
  splits <- unique(sce$treatment)
  names(splits) <- splits
  pb_list <- splits |>
    lapply(function(tr) sce[, sce$treatment == tr]) |>
    lapply(function(sub) {
      pb <- rhapsodykit::calculate_pseudo_bulk(
        sub, type, fun = fun, by = c("time", "cellType"))
      as.list(pb@assays@data)
    })

  pb_list <- unlist(pb_list, recursive = FALSE, use.names = TRUE) |>
    purrr::imap(function(mtx, name) {
      colnames(mtx) <- paste(name, colnames(mtx), sep = ".")
      mtx
    })

  do.call(cbind, pb_list)
}

treatment_colors <- c(
  "MOCK" = "blue",
  "PNR2" = "red",
  "TR4"  = "forestgreen"
)

formatCosgTable <- function(cosg_res) {
  lapply(names(cosg_res$names), function(cls) {
    data.frame(
      cluster = cls,
      names = cosg_res$names[[cls]],
      scores = cosg_res$scores[[cls]]
    )
  }) |> dplyr::bind_rows()
}

ggpreview <- function(ggobj, ext = ".png",...) {
  temp.file <- paste0(tempfile(), ext)
  ggplot2::ggsave(filename = temp.file, plot = ggobj, ...)
  utils::browseURL(temp.file)
}

pltpreview <- function(
    fig, ext = ".png",
    width = 7, height = 7, dpi = 300, ...) {
  temp.file <- paste0(tempfile(), ext)
  fig$set_size_inches(width, height)
  fig$savefig(temp.file, dpi = dpi, ...)
  utils::browseURL(temp.file)
}

plotpreview <- function(x, width = 7, height = 7, res = 300, ...) {
  temp.file <- paste0(tempfile(), ".png")
  png(filename = temp.file, units = "in",
      width = width, height = height, res = res, ...)
  print(x)
  dev.off()
  utils::browseURL(temp.file)
  invisible()
}

rotate_x_labels <- function(angle = 45, hjust = 1, ...) {
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = angle, hjust = hjust, ...))
}

verticalize_x_labels <- function(angle = 90, hjust = 1, vjust = 0.5, ...) {
  rotate_x_labels(angle = angle, hjust = hjust, vjust = vjust, ...)
}

center_plot_title <- function(...) {
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, ...))
}

remove_axis <- function(...) {
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    ...)
}

remove_x_axis <- function(...) {
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_blank(),
    ...)
}

remove_legend <- function(...) {
  ggplot2::theme(
    legend.position = "none", validate = TRUE, ...)
}

remove_legend_title <- function(...) {
  ggplot2::theme(
    legend.title = ggplot2::element_blank(), ...)
}

empty_strip <- function(...) {
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    ...
  )
}

remove_strip <- function(...) {
  ggplot2::theme(
    strip.text = ggplot2::element_blank() , 
    strip.background = ggplot2::element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0) , units = "lines"),
    ...)
}

legend_override <- function(key, values, ...) {
  args <- list()
  args[[key]] <- ggplot2::guide_legend(
    override.aes = values, ...)
  do.call(ggplot2::guides, args)
}

font_plot_tag <- function(size = 11, ...) {
  ggplot2::theme(plot.tag = ggplot2::element_text(
    size = size, face = "bold", ...))
}

font_text_title <- function(...) {
  ggplot2::theme(text = ggplot2::element_text(...),
                 title = ggplot2::element_text(...))
}

theme_dimred <- function(
    xlength = 0.3, ylength = 0.3, 
    arrow = grid::arrow(
      angle = 15, length = grid::unit(0.15, "inches"),
      type = "closed"),
    ...) {
  ggplot2::theme(
    axis.line = ggplot2::element_blank(), 
    axis.ticks = ggplot2::element_blank(), 
    axis.text = ggplot2::element_blank(),
    axis.line.x.bottom = tidydr::element_line2(
      id = 1, xlength = xlength, arrow = arrow),
    axis.line.y.left = tidydr::element_line2(
      id = 2, ylength = ylength, arrow = arrow), 
    axis.title = ggplot2::element_text(hjust=0.1),
    ...
  )
}

download_table <- function(df, table_id = "table", ...) {
  htmltools::tagList(
    htmltools::tags$button(
      "Download as CSV",
      onclick = sprintf(
        "Reactable.downloadDataCSV('%s', '%s.csv')",
        table_id, table_id)
    ),
    reactable::reactable(df, elementId = table_id, ...)
  )
}


