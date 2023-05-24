# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

#' Run Clustering Pipeline Used by This Project
RunSeuratPipeline <- function(
    seurat_object, verbose = FALSE, nfeatures = 2000,
    num.pc = 20, resolution = 0.8, k.param = 30) {

  seurat_object |>
    Seurat::NormalizeData(verbose = verbose) |>
    Seurat::FindVariableFeatures(
      selection.method = "vst", nfeatures = nfeatures, verbose = verbose) |>
    Seurat::ScaleData(verbose = verbose) |>
    Seurat::RunPCA(features = NULL, npcs = num.pc, verbose = verbose) |>
    harmony::RunHarmony(
      group.by.vars = c("time", "sample"), theta = c(4, 2),
      kmeans_init_nstart = 20, kmeans_init_iter_max = 100,
      epsilon.cluster = -Inf, epsilon.harmony = -Inf,
      plot_convergence = FALSE, verbose = verbose) |>
    Seurat::FindNeighbors(
      reduction = "harmony", dims = seq_len(num.pc),
      k.param = k.param, verbose = verbose, force.recalc = TRUE) |>
    Seurat::FindClusters(resolution = resolution, verbose = verbose)
}

RandomSubsetData<- function(object, rate, random.subset.seed = NULL, ...){
  ncells <- nrow(object@meta.data)
  ncells.subsample <- round(ncells * rate)
  set.seed(random.subset.seed)
  selected.cells <- sample(colnames(object), ncells.subsample)
  subset(object, cells =  selected.cells, ...)
}
