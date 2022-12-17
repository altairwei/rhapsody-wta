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
          searchable = TRUE,
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
      searchable = TRUE,
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
      searchable = TRUE,
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
    method = c("sum", "gsva", "moduleScore"), ...) {
  
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
     method = "gsva", ...),
    moduleScore = calculateModuleScore(
      expr, pathways, ...)
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

lengthCheck <- function(values, cutoff = 0) {
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}

caseMatch <- function(search, match) {
  search.match <- sapply(
    X = search,
    FUN = function(s) {
      return(grep(
        pattern = paste0('^', s, '$'),
        x = match,
        ignore.case = TRUE,
        perl = TRUE,
        value = TRUE
      ))
    }
  )
  return(unlist(x = search.match))
}

#' Calculate module scores for feature expression programs in single cells
#'
#' Inspired by https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/
#'
#' Calculate the average expression levels of each program (cluster) on single
#' cell level, subtracted by the aggregated expression of control feature sets.
#' All analyzed features are binned based on averaged expression, and the
#' control features are randomly selected from each bin.
#'
#' @param object A genes × cells expression matrix
#' @param features A (named) list of vectors of features for expression programs;
#' each entry should be a vector of feature names
#' @param pool List of features to check expression levels against, defaults to
#' \code{rownames(x = object)}
#' @param nbin Number of bins of aggregate expression levels for all
#' analyzed features
#' @param ctrl Number of control features selected from the same bin per
#' analyzed feature
#' @param name Name for the expression programs; If \code{features} is not named,
#' \{name} will append a number to the end for each entry in \code{features} (eg.
#' if \code{features} has three programs, the results will be stored as \code{name1},
#' \code{name2}, \code{name3}, respectively).
#' @param seed Set a random seed. If NULL, seed is not set.
#'
#' @return Returns a modules × cells matrix of module scores
#'
#' @references Tirosh et al, Science (2016)
calculateModuleScore <- function(
    object,
    features,
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    assay = NULL,
    name = 'Cluster',
    seed = 1,
    ...
) {
  if (!is.null(seed))
    set.seed(seed)

  if (is.null(x = features)) {
    stop("Missing input feature list")
  }

  features <- lapply(
    X = features,
    FUN = function(x) {
      missing.features <- setdiff(x, y = rownames(object))
      if (length(missing.features) > 0) {
        warning(
          "The following features are not present in the object: ",
          paste(missing.features, collapse = ", "),
          call. = FALSE,
          immediate. = TRUE
        )
      }
      return(intersect(x, rownames(object)))
    }
  )

  cluster.length <- length(features)

  if (!all(lengthCheck(features))) {
    warning(paste(
      'Could not find enough features in the object from the following feature lists:',
      paste(names(which(!lengthCheck(features)))),
      'Attempting to match case...'
    ))
    features <- lapply(
      X = features.old,
      FUN = caseMatch,
      match = rownames(object)
    )
  }

  if (!all(lengthCheck(features))) {
    stop(paste(
      'The following feature lists do not have enough features present in the object:',
      paste(names(x = which(x = !lengthCheck(features)))),
      'exiting...'
    ))
  }

  if (is.null(pool))
    pool = rownames(object)

  # Get the average expression across all cells (named vector)
  data.avg <- Matrix::rowMeans(object[pool, ])
  # Order genes from lowest average expression to highest average expression
  data.avg <- data.avg[order(data.avg)]

  # Use ggplot2's cut_number function to make n groups with (approximately) equal
  # numbers of observations. The 'rnorm(n = length(data.avg))/1e+30' part adds a
  # tiny bit of noise to the data, presumably to break ties.
  data.cut <- ggplot2::cut_number(
    x = data.avg + rnorm(length(data.avg)) / 1e+30,
    n = nbin, labels = FALSE, right = FALSE)

  # Set the names of the cuts as the gene names
  names(data.cut) <- names(data.avg)

  # Create an empty list the same length as the number of input gene sets. This will
  # contain the names of the control genes
  ctrl.use <- vector(mode = "list", length = cluster.length)

  for (i in 1:cluster.length) {
    # Get the gene names from the input gene set as a character vector  
    features.use <- features[[i]]

    # Loop through the provided genes (1:num_genes) and for each gene, find ctrl
    # (default=100) genes from the same expression bin (by looking in data.cut):
    for (j in 1:length(features.use)) {
      # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin
      # number. We then sample `ctrl` genes from that bin without replacement and
      # add the gene names to ctrl.use.
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(sample(
          data.cut[which(data.cut == data.cut[features.use[j]])],
          size = ctrl, replace = FALSE)
        )
      )
    }
  }

  # Remove any repeated gene names - even though we set replace=FALSE when we
  # sampled genes from the same expression bin, there may be more than two genes
  # in our input gene list that fall in the same expression bin, so we can end
  # up sampling the same gene more than once.
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)

  ## Get control gene scores

  # Create an empty matrix with dimensions;
  # number of rows equal to the number of gene sets (just one here)
  # number of columns equal to number of cells in input Seurat object
  ctrl.scores <- matrix(data = numeric(length = 1L),
                        nrow = length(ctrl.use),
                        ncol = ncol(object))

  # Loop through each provided gene set and add to the empty matrix the mean
  # expression of the control genes in each cell
  for (i in 1:length(ctrl.use)) {
    # Get control gene names as a vector  
    features.use <- ctrl.use[[i]]
    # For each cell, calculate the mean expression of *all* of the control genes 
    ctrl.scores[i, ] <- Matrix::colMeans(object[features.use,])
  }

  ## Get scores for input gene sets

  # Similar to the above, create an empty matrix
  features.scores <- matrix(data = numeric(length = 1L),
                            nrow = cluster.length,
                            ncol = ncol(x = object))

  # Loop through input gene sets and calculate the mean expression of these
  # genes for each cell
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- object[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }

  # Subtract the control scores from the feature scores - the idea is that
  # if there is no enrichment of the genes in the geneset in a cell, then the
  # result of this subtraction should be ~ 0
  features.scores.use <- features.scores - ctrl.scores

  # Name the result the "name" variable + whatever the position the geneset was
  # in the input list, e.g. "Cluster1"
  module_names <- names(features)
  if (!is.null(module_names) && length(module_names) == cluster.length)
    rownames(features.scores.use) <- names(features)
  else
    rownames(features.scores.use) <- paste0(name, 1:cluster.length)

  # Give the rows of the matrix, the names of the cells
  colnames(features.scores.use) <- colnames(object)

  features.scores.use
}

#' Remove duplicated GO term
enrichDistinct <- function(x) {
  # Split core_enrichment or genes into list, then sort them by alphabet.
  if (inherits(x, "gseaResult")) {
    dups <- dplyr::mutate(x@result, geneList = lapply(
        strsplit(core_enrichment, "/"), sort)) |>
      dplyr::group_by(geneList) |>
      tidyr::nest(data = !geneList) |>
      dplyr::mutate(n = sapply(data, nrow)) |>
      dplyr::filter(n > 1) |>
      dplyr::mutate(dups = lapply(data, function(df) {
        df$IDNUM <- as.integer(substring(df$ID, 4))
        df$ID[df$IDNUM != max(df$IDNUM)]
      })) |>
      dplyr::pull(dups) |>
      unlist()

    dplyr::filter(x, !ID %in% dups)
  } else {
    x
  }
}
