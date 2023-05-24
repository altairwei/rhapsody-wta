source("renv/activate.R")

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

output_file <- snakemake@output[[1]]
objfiles <- snakemake@input

obj_list <- lapply(snakemake@input, readRDS)
names(obj_list) <- sapply(obj_list, Seurat::Project)

if (is.null(names(obj_list)))
  stop("Seurat object must have a project name.")

obj_merged <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list)
)

saveRDS(obj_merged, file = output_file)
