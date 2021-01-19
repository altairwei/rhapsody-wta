save_plot <- function(filename, plot, width = 7, height = 7, ...) {
  message(sprintf("Saving %d x %d in image: %s", width, height, filename))
  ggplot2::ggsave(filename, plot, width = width, height = height, ...)
}

save_png <- function(filename, plot, width = 7, height = 7, dpi = 300) {
  message(sprintf("Saving %d x %d in image: %s", width, height, filename))
  png(filename, width = width, height = height, units = "cm", res = dpi)
  print(plot)
  dev.off()
}

read_raw_csv <- function(base_dir = ".") {
  if (!dir.exists(base_dir)) {
    stop("Directory provided does not exist")
  }

  matrix_loc <- Sys.glob(file.path(base_dir, "*_RSEC_MolsPerCell.csv"))

  if (length(matrix_loc) != 1) {
    stop("`*_RSEC_MolsPerCell.csv` missing or more than one file was found.")
  }

  message(sprintf("Reading %s", matrix_loc))
  df <- readr::read_csv(matrix_loc, comment = "#", progress = TRUE,
    col_types = readr::cols(
      Cell_Index = readr::col_character()
    )
  )

  message("Matrix constructing")
  m <- data.matrix(df[, -1])
  rownames(m) <- df[["Cell_Index"]]
  m <- t(m)
  mtx <- Matrix::Matrix(m)

  rm(df, m)
  gc()

  mtx
}

read_mtx <- function(base_dir = ".") {
  if (!dir.exists(base_dir)) {
    stop("Directory provided does not exist")
  }

  matrix_loc <- Sys.glob(
    file.path(base_dir, "*_Expression_Matrix.mtx"))

  if (length(matrix_loc) != 1) {
    stop("`*_Expression_Matrix.mtx` missing or more than one file was found.")
  }

  colnames_loc <- sprintf("%s.colnames", matrix_loc)
  rownames_loc <- sprintf("%s.rownames", matrix_loc)
  stopifnot(file.exists(colnames_loc, rownames_loc))

  message(sprintf("Reading %s", matrix_loc))
  m <- Matrix::readMM(matrix_loc)
  colnames(m) <- readr::read_lines(colnames_loc)
  rownames(m) <- readr::read_lines(rownames_loc)

  m
}

read_rhapsody_wta <- function(base_dir, use_mtx = FALSE) {
  expr_matrix <- NULL

  if (isTRUE(use_mtx)) {
    expr_matrix <- read_mtx(base_dir)
  } else {
    expr_matrix <- read_raw_csv(base_dir)
  }

  expr_matrix
}