source("renv/activate.R")
snakemake@source("common.R")

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

seurat_obj <- readRDS(snakemake@input[[1]])
output_obj <- snakemake@output[["obj"]]
output_df <- snakemake@output[["df"]]

k<- as.numeric(snakemake@wildcards[["k"]])
resolution <- as.numeric(snakemake@wildcards[["resolution"]])
npc <- as.numeric(snakemake@wildcards[["pc"]])

verbose <- snakemake@params[["verbose"]] %||% FALSE
nfeatures <- snakemake@params[["nfeatures"]] %||% 2000L

seurat_obj <- RunSeuratPipeline(
  seurat_obj, k.param = k, num.pc = npc, resolution = resolution,
  verbose = verbose, nfeatures = nfeatures)

res_df <- tibble::tibble(
  pc = npc,
  resolution = resolution,
  k_param = k,
  original_ident_full = list(Seurat::Idents(seurat_obj))
)

saveRDS(seurat_obj, file = output_obj)
saveRDS(res_df, file = output_df)
