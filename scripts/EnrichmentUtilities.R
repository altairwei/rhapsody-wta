library(magrittr)
suppressMessages(library(ggplot2))


#' Perform ORA enrichment analysis for DS results
#'
#' @param ds DS results
#' @param contrasts Which contrast to pull
#' @param clusters Which clusters to compare
#' @param ont Gene ontology to display
#' @param simplify Simplify GO term
#'
#' @return A list of two elements: \code{enr} and \code{fc}
#'
performDSEnrichORA <- function(
  ds, contrasts, clusters,
  ont = c("BP", "MF", "CC", "ALL"),
  simplify = TRUE
) {
  ont <- match.arg(ont)
  
  # We use all detected genes in DS analysis of given cluster as background
  background <- rhapsodykit::diff_state_pull(
    ds, contrasts, clusters, "gene")
  
  sig_genes <- ds %>%
    rhapsodykit::diff_state_significant()
  genes <- sig_genes %>%
    rhapsodykit::diff_state_pull(contrasts, clusters, "gene")
  fc <- sig_genes %>%
    rhapsodykit::diff_state_pull(contrasts, clusters, "logFC")
  names(fc) <- genes
  
  enr <- clusterProfiler::enrichGO(
    genes,
    org.Taestivum.iwgsc.db,
    keyType = "GID",
    ont = ont,
    universe = background
  )
  
  if (simplify && ont != "ALL")
    enr <- clusterProfiler::simplify(enr, cutoff = 0.5)
  
  list(
    enr = enr,
    fc = fc
  )
}

#' Compare ORA enrichment results
#'
#' @inheritParams performDSEnrichORA
#'
#' @return A list of two elements: \code{enr} and \code{fc}
#'
performDSCompareORA <- function(
  ds, contrasts, clusters,
  ont = c("BP", "MF", "CC", "ALL"),
  simplify = TRUE
) {
  ont <- match.arg(ont)
  
  if (is.null(names(contrasts)))
    names(contrasts) <- contrasts
  if (is.null(names(clusters)))
    names(clusters) <- clusters
  
  background <- purrr::imap(clusters, function(clr, clr_name) {
    purrr::imap(contrasts, function(con, con_name) {
      ds %>% rhapsodykit::diff_state_pull(con, clr, "gene")
    })
  }) %>%
    unlist() %>%
    unique()
  
  compare_df <- purrr::imap(clusters, function(clr, clr_name) {
    purrr::imap(contrasts, function(con, con_name) {
      df <- ds %>%
        rhapsodykit::diff_state_significant() %>%
        rhapsodykit::diff_state_pull(con, clr, c("gene", "logFC"))
      df$gene
      data.frame(
        gene = df$gene,
        contrast = con_name,
        cluster = clr_name
      )
    })
  }) %>%
    unlist(recursive = FALSE) %>%
    dplyr::bind_rows()
  
  compare_df$cluster <- factor(compare_df$cluster, levels = names(clusters))
  compare_df$contrast <- factor(compare_df$contrast, levels = names(contrasts))
  
  comp_enr <- clusterProfiler::compareCluster(
    gene ~ cluster + contrast,
    data = compare_df,
    fun = "enrichGO",
    OrgDb = org.Taestivum.iwgsc.db,
    keyType = "GID",
    ont = ont,
    universe = background
  )
  
  if (simplify && ont != "ALL")
    comp_enr <- clusterProfiler::simplify(comp_enr, cutoff = 0.5)
  
  comp_enr
}

performDSEnrichGSEA <- function(
  ds, contrasts, clusters,
  ont = c("BP", "MF", "CC", "ALL"),
  simplify = TRUE
) {
  ont <- match.arg(ont)
  
  ## 我们不需要依据 FDR 来筛选基因，只需要 logFC 来排序
  gsea_df <- ds %>%
    rhapsodykit::diff_state_pull(contrasts, clusters, c("gene", "logFC"))
  
  ## feature 1: numeric vector
  geneList <- gsea_df[, 2]
  
  ## feature 2: named vector
  names(geneList) <- as.character(gsea_df[, 1])
  
  ## feature 3: decreasing order
  geneList <- sort(geneList, decreasing = TRUE)
  
  gsea <- clusterProfiler::gseGO(
    geneList,
    ont = ont,
    OrgDb = org.Taestivum.iwgsc.db,
    keyType = "GID"
  )
  
  if (simplify && ont != "ALL")
    gsea <- clusterProfiler::simplify(gsea)
  
  list(
    gsea = gsea,
    fc = geneList
  )
  
}

performDSCompareGSEA <- function(
  ds, contrasts, clusters,
  ont = c("BP", "MF", "CC", "ALL"),
  simplify = TRUE
) {
  ont <- match.arg(ont)
  
  if (is.null(names(contrasts)))
    names(contrasts) <- contrasts
  if (is.null(names(clusters)))
    names(clusters) <- clusters
  
  compare_df <- purrr::imap(clusters, function(clr, clr_name) {
    purrr::imap(contrasts, function(con, con_name) {
      df <- ds %>%
        rhapsodykit::diff_state_pull(con, clr, c("gene", "logFC"))
      geneList <- df[, 2]
      names(geneList) <- as.character(df[, 1])
      geneList <- sort(geneList, decreasing = TRUE)
      data.frame(
        gene = names(geneList),
        fc = geneList,
        cluster = clr_name,
        contrast = con_name
      )
    })
  }) %>%
    unlist(recursive = FALSE) %>%
    dplyr::bind_rows()
  
  compare_df$cluster <- factor(compare_df$cluster, levels = names(clusters))
  compare_df$contrast <- factor(compare_df$contrast, levels = names(contrasts))
  
  comp_enr <- clusterProfiler::compareCluster(
    gene | fc ~ cluster + contrast,
    data = compare_df,
    OrgDb = org.Taestivum.iwgsc.db,
    fun = "gseGO",
    keyType = "GID",
    ont = ont
  )
  
  if (simplify && ont != "ALL")
    comp_enr <- clusterProfiler::simplify(comp_enr, cutoff = 0.5)
  
  comp_enr
}

#' Perform comparative ORA for cell type marker genes
#'
#' @inheritParams performDSEnrichORA
#' @param markerlist A list of character vectors.
#' @param clusters A (named) character vector which represents
#' selected cell types.
#' @param background Background genes used for ORA.
#' @return Results of \code{clusterProfiler::compareCluster}
#'
performMarkerCompareORA <- function(
    markerlist, clusters, background,
    ont = c("BP", "MF", "CC", "ALL"),
    simplify = TRUE
) {
  ont <- match.arg(ont)
  
  if (is.null(names(clusters)))
    names(clusters) <- clusters
  
  compare_df <- purrr::imap(clusters, function(clr, clr_name) {
    df <- data.frame(
      gene = markerlist[[clr]],
      cluster = rep(clr_name, length(markerlist[[clr]]))
    )
    df
  }) %>%
    dplyr::bind_rows()
  
  compare_df$cluster <- factor(compare_df$cluster, levels = names(clusters))
  
  comp_enr <- clusterProfiler::compareCluster(
    gene ~ cluster,
    data = compare_df,
    fun = "enrichGO",
    OrgDb = org.Taestivum.iwgsc.db,
    keyType = "GID",
    ont = ont,
    universe = background
  )
  
  if (simplify && ont != "ALL")
    comp_enr <- clusterProfiler::simplify(comp_enr, cutoff = 0.5)
  
  comp_enr
}

#' Print Compare ORA Enrichment Results by Cluster
#'
#' @param enr Comparative ORA enrichment results.
#' @param heading Markdown heading symbols.
printCompareORAByCluster <- function(enr, heading) {
  df <- enr@compareClusterResult |>
    dplyr::select(
      cluster, ID, Description,
      Count, p.adjust, geneID) |>
    dplyr::mutate(
      p.adjust = format(p.adjust, scientific = TRUE, digits = 2),
      geneList = sapply(openssl::md5(geneID), function(x) substr(x, 1, 9))
    ) |>
    dplyr::group_by(cluster)
  
  df_list <- dplyr::group_split(df)
  names(df_list) <- dplyr::group_keys(df)$cluster
  
  for (clr in names(df_list)) {
    cat(heading, clr, "\n\n")
    print(
      htmltools::tagList(
        htmltools::tags$button(
          "Download as CSV",
          onclick = sprintf(
            "Reactable.downloadDataCSV('%s', '%s')",
            clr, clr)
        ),
        reactable::reactable(
          df_list[[clr]],
          elementId = clr,
          columns = list(
            "cluster" = reactable::colDef(maxWidth = 80),
            "ID" = reactable::colDef(maxWidth = 120),
            "Description" = reactable::colDef(minWidth = 100),
            "Count" = reactable::colDef(maxWidth = 80),
            "p.adjust" = reactable::colDef(
              maxWidth = 100, format = reactable::colFormat(digits = 2)),
            "geneID" = reactable::colDef(minWidth = 250, show = FALSE),
            "geneList" = reactable::colDef(
              maxWidth = 140,
              details = function(row_idx) {
                reactable::reactable(
                  data.frame(geneID = df_list[[clr]][row_idx, "geneID"]),
                  outlined = FALSE, fullWidth = TRUE)
              }
            )
          )
        )
      )
    )
    cat("\n\n")
  }
}

performDACompareORA <- function(da, sets,
                            ont = c("BP", "MF", "CC", "ALL"),
                            simplify = TRUE,
                            logFC = c("both", "up", "down")
) {
  stopifnot(is.list(sets) && !is.null(names(sets)) && length(sets) > 1)
  logFC = match.arg(logFC)

  compare_list <- sets %>%
    lapply(function(set) da[[ set[[1]] ]][[ set[[2]] ]])
  
  comp_enr <- clusterProfiler::compareCluster(
    compare_list, fun = "enrichGO",
    OrgDb = org.Taestivum.iwgsc.db,
    keyType = "GID",
    ont = ont)

  if (simplify && ont != "ALL")
    comp_enr <- clusterProfiler::simplify(comp_enr)

  comp_enr
}

#' Print GSEA Results with Reactable
#'
#' @param gsea GSEA result object
#' @return Reactable object
printGSEATable <- function(gsea) {
  enrichdf <- gsea@result |>
    tibble::remove_rownames() |>
    dplyr::mutate(
      p.adjust = format(p.adjust, scientific = TRUE, digits = 2),
      Count = sapply(strsplit(core_enrichment, "/"), length),
      geneList = sapply(openssl::md5(core_enrichment), function(x) substr(x, 1, 9))) |>
    dplyr::select(ID, Description, NES,
                  Count, p.adjust, geneList, core_enrichment)

  htmltools::tagList(
    htmltools::tags$button(
      "Download as CSV",
      onclick = "Reactable.downloadDataCSV('GSEATable', 'GSEATable.csv')"
    ),
    reactable::reactable(
      enrichdf,
      elementId = "GSEATable",
      columns = list(
        "ID" = reactable::colDef(maxWidth = 120),
        "Description" = reactable::colDef(minWidth = 100),
        "NES" = reactable::colDef(maxWidth = 80,
                                  format = reactable::colFormat(digits = 2)),
        "Count" = reactable::colDef(maxWidth = 80),
        "p.adjust" = reactable::colDef(maxWidth = 100),
        "core_enrichment" = reactable::colDef(show = FALSE),
        "geneList" = reactable::colDef(
          maxWidth = 140,
          details = function(row_idx) {
            reactable::reactable(
              data.frame(geneID = enrichdf[row_idx, "core_enrichment"]),
              outlined = FALSE, fullWidth = TRUE)
          }
        )
      )
    )
  )
}

printCompareORATable <- function(ora) {
  enrichdf <- ora@compareClusterResult |>
    tibble::remove_rownames() |>
    dplyr::mutate(
      p.adjust = format(p.adjust, scientific = TRUE, digits = 2),
      geneList = sapply(openssl::md5(geneID), function(x) substr(x, 1, 9))) |>
    dplyr::select(Cluster, ID, Description, GeneRatio,
                  Count, p.adjust, geneID, geneList)

  htmltools::tagList(
    htmltools::tags$button(
      "Download as CSV",
      onclick = "Reactable.downloadDataCSV('ORATable', 'ORATable.csv')"
    ),
    reactable::reactable(
      enrichdf,
      elementId = "ORATable",
      columns = list(
        "Cluster" = reactable::colDef(maxWidth = 120),
        "ID" = reactable::colDef(maxWidth = 120),
        "Description" = reactable::colDef(minWidth = 100),
        "GeneRatio" = reactable::colDef(maxWidth = 100),
        "Count" = reactable::colDef(maxWidth = 80),
        "p.adjust" = reactable::colDef(
          maxWidth = 100, format = reactable::colFormat(digits = 2)),
        "geneID" = reactable::colDef(show = FALSE),
        "geneList" = reactable::colDef(
          maxWidth = 140,
          details = function(row_idx) {
            reactable::reactable(
              data.frame(geneID = enrichdf[row_idx, "geneID"]),
              outlined = FALSE, fullWidth = TRUE)
          }
        )
      )
    )
  )
}


#' Calculate Pathway Activity
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param pathways Gene sets provided as a \code{list} object.
#' @param annotation A \code{data.frame} object of pathway annotation.
#' Rownames should be the pathway id same to names of \code{pathways} list.
#' @param by_exprs_values Which kind of expression used for calculation.
#' @param method The way to calculate activity. \code{sum} means summed
#' expression of genes in each pathway. \code{gsva} means GSVA enrichment
#' scores.
#' @param ... Arguments passed to \code{\link[GSVA]{gsva}}
#' @return A \code{SingleCellExperiment} object
calcPathwayActivity <- function(
    sce, pathways, annotation = NULL,
    by_exprs_values = "logcounts",
    method = c("sum", "gsva"), ...) {
  
  expr <- assay(sce, by_exprs_values)
  
  method <- match.arg(method)
  pathway_expr <- switch(method,
    sum = t(vapply(
     X = pathways,
     FUN = function(genes)
       Matrix::colSums(expr[genes,]),
     FUN.VALUE = numeric(ncol(expr))
    )),
    gsva = GSVA::gsva(
     expr[unlist(pathways), ], pathways,
     method = "gsva", ...)
  )
  
  SingleCellExperiment(
    assays = list(activity = pathway_expr),
    reducedDims = reducedDims(sce),
    rowData = annotation[rownames(pathway_expr),]
  )
}

#' Plot Pathway Activation on Cell Embeddings
#'
#' @param sce SingleCellExperiment object with pathway ID as rows.
#' @param id Pathway ID to plot.
#' @param order Order cells according to activity score.
#' @param ... Arguments passed to \code{\link[scater]{plotReducedDim}}
#'
#' @return A ggplot object
plotPathwayActivation <- function(sce, id, order = TRUE, ...) {
  #Require scater >= 1.26.0
  scater::plotReducedDim(
    sce,
    by_exprs_values = "activity",
    colour_by = id, order_by = id, ...) +
    ggplot2::scale_color_gradient(low = "grey", high = "blue") +
    ggplot2::ggtitle(rowData(sce)[id, "Description"])
}
