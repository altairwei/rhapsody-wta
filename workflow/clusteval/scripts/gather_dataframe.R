source("renv/activate.R")
snakemake@source("common.R")

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

rdss <- snakemake@input[["rds"]]

dat.list <- lapply(rdss, readRDS)
gather_idents <- do.call(dplyr::bind_rows, dat.list)
saveRDS(gather_idents, file = snakemake@output[[1]])
