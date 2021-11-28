#!/usr/bin/env Rscript

workflowr::wflow_build(
    c(
        "analysis/sample_integration.Rmd",
        "analysis/cell_annotation.Rmd",
        "analysis/differential_state_analysis.Rmd",
        "analysis/about.Rmd",
        "analysis/license.Rmd",
        "analysis/index.Rmd"
    ),
    make = FALSE
)