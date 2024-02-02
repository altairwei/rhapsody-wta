source("scripts/LoadUtils.R", chdir = TRUE)
sce_os_all <- readRDS(Sys.glob("results/ObjectCache/TrajectoryInference/sce_all_slingshot_os_all_*.rds"))
patternRes_os_all <- readRDS(Sys.glob("results/ObjectCache/TrajectoryInference/tradeseq_patternRes_os_all_*.rds"))

patternGenes <- patternRes_os_all |>
  tibble::rownames_to_column("gene") |>
  dplyr::mutate(fdr = p.adjust(pvalue, "fdr")) |>
  dplyr::filter(fdr <= 0.01) |>
  dplyr::pull(gene) |>
  unique()

tradeSce_os_cond <- xfun::cache_rds(
  file = "tradeseq_gam_os_condition.rds",
  dir = "results/ObjectCache/TrajectoryInference/",
  expr = tradeSeq::fitGAM(
    # fitGAM needs counts rather than logcounts, see
    # the source code of SCE-version fitGAM.
    counts = sce_os_all,
    conditions = factor(sce_os_all$treatment),
    genes = patternGenes,
    nknots = 8,
    sce = TRUE,
    verbose = TRUE,
    parallel = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(workers = 1)
  )
)
