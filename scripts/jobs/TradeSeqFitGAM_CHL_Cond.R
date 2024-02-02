source("scripts/LoadUtils.R", chdir = TRUE)
sce_chl_all <- readRDS(Sys.glob("results/ObjectCache/TrajectoryInference/sce_all_slingshot_chl_all_*.rds"))
patternRes_chl_all <- readRDS(Sys.glob("results/ObjectCache/TrajectoryInference/tradeseq_patternRes_chl_all_*.rds"))

patternGenes <- patternRes_chl_all |>
  tibble::rownames_to_column("gene") |>
  dplyr::mutate(fdr = p.adjust(pvalue, "fdr")) |>
  dplyr::filter(fdr <= 0.01) |>
  dplyr::pull(gene) |>
  unique()

tradeSce_chl_cond <- xfun::cache_rds(
  file = "tradeseq_gam_chl_condition.rds",
  dir = "results/ObjectCache/TrajectoryInference/",
  expr = tradeSeq::fitGAM(
    # fitGAM needs counts rather than logcounts, see
    # the source code of SCE-version fitGAM.
    counts = sce_chl_all,
    conditions = factor(sce_chl_all$treatment),
    genes = patternGenes,
    nknots = 8,
    sce = TRUE,
    verbose = TRUE,
    parallel = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(workers = 1)
  )
)
