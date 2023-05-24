source("renv/activate.R")
snakemake@source("common.R")

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

seurat_obj <- readRDS(snakemake@input[[1]])
output_file <- snakemake@output[[1]]

k <- as.numeric(snakemake@wildcards[["k"]])
resolution <- as.numeric(snakemake@wildcards[["resolution"]])
npc <- as.numeric(snakemake@wildcards[["pc"]])
run_id <- snakemake@wildcards[["run_id"]]

verbose <- snakemake@params[["verbose"]] %||% FALSE
nfeatures <- snakemake@params[["nfeatures"]] %||% 2000L

subset_seurat_obj <- RandomSubsetData(seurat_obj, rate = snakemake@params[["rate"]])
original_ident <- Seurat::Idents(subset_seurat_obj)

## after reprocessing, the ident slot will be updated with the new cluster id
subset_seurat_obj <- RunSeuratPipeline(
  subset_seurat_obj, k.param = k, num.pc = npc, resolution = resolution,
  verbose = verbose, nfeatures = nfeatures)

res <- tibble::tibble(
  pc = npc, resolution = resolution, k_param = k,
  original_ident = list(original_ident),
  recluster_ident = list(Seurat::Idents(subset_seurat_obj)),
  round = run_id)

saveRDS(res, file = output_file)

## make sure it is not empty file
info <- file.info(output_file)
if (info$size == 0) {
  quit(status = 1)
}
