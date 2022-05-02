#!/usr/bin/env Rscript

library(magrittr)
library(optparse)

parser <- OptionParser()
arguments <- parse_args2(parser)

clean_images <- function(htmlfile) {
  message("Checking ", fs::path_file(htmlfile))

  htmldir <- fs::path_dir(htmlfile)

  images <- rvest::read_html(htmlfile) %>%
    rvest::html_elements("img") %>%
    rvest::html_attr("src")

  if (length(images) == 0) {
    message("No image")
    return()
  }

  images <- paste(htmldir, images, sep = "/")

  figdir <- fs::path_join(c(
    fs::path_dir(htmlfile),
    "figure",
    fs::path_file(fs::path_ext_set(htmlfile, "Rmd"))
  ))

  img_to_delete <- setdiff(fs::dir_ls(figdir), images)

  if (length(images) == length(img_to_delete)) {
    message("Not to delete all images, skipping.")
    return()
  }

  if (length(img_to_delete) == 0) {
    message("No image to delete.")
  } else {
    message("Deleting ", length(img_to_delete), " images.")
  }

  fs::file_delete(img_to_delete)
}

for (file in arguments$args) {
  clean_images(file)
}