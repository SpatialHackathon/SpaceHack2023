#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created script

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(
    c("-o", "--out_dir"),
    type = "character", default = NULL,
    help = "Output directory to write files to."
  )
)

description <- "Load data for LIBD DLPFC (http://research.libd.org/spatialLIBD/)."

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

out_dir <- opt$out_dir

## Your code goes here
suppressPackageStartupMessages(library(spatialLIBD))
suppressPackageStartupMessages(library(magrittr))

write_SpatialExperiment_to_folder <- function(
    spe, path, obs_col, label_col = "label", assay_name = "counts") {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  colData(spe)[label_col] %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(!!as.symbol(label_col))) %>%
    write.table(file.path(path, "labels.tsv"), sep = "\t", col.names = NA, quote = FALSE)

  colData(spe)[obs_col] %>%
    as.data.frame() %>%
    write.table(file.path(path, "observations.tsv"), sep = "\t", col.names = NA, quote = FALSE)

  rowData(spe) %>%
    as.data.frame() %>%
    write.table(file.path(path, "features.tsv"), sep = "\t", col.names = NA, quote = FALSE)

  coords <- spatialCoords(spe)
  mode(coords) <- "integer"
  as.data.frame(coords) %>%
    dplyr::rename(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres") %>%
    write.table(file.path(path, "coordinates.tsv"), sep = "\t", col.names = NA, quote = FALSE)

  assay(spe, assay_name) %>%
    t() %>%
    Matrix::writeMM(file.path(path, "counts.mtx"))
}

spe <- fetch_data("spe")

keep_cols <- c("sample_id", "subject", "position", "replicate", "discard", "spatialLIBD", "array_row", "array_col")

colData(spe) <- colData(spe)[, keep_cols]
colnames(colData(spe))[colnames(colData(spe)) == "array_row"] <- "row"
colnames(colData(spe))[colnames(colData(spe)) == "array_col"] <- "col"
colnames(colData(spe))[colnames(colData(spe)) == "spatialLIBD"] <- "label"

keep_rows <- c("gene_version", "gene_name", "source", "gene_biotype")
rowData(spe) <- rowData(spe)[, keep_rows]

patients <- unique(colData(spe)$subject)
for (patient in patients) {
  patient_spe <- spe[, spe$subject == patient]
  samples <- unique(colData(patient_spe)$sample_id)
  for (sample in samples) {
    spe_sample <- patient_spe[, patient_spe$sample_id == sample]
    colData(spe_sample) <- colData(spe_sample)[, c("label", "row", "col")] # "discard"
    write_SpatialExperiment_to_folder(
      spe_sample,
      file.path(out_dir, paste(patient, sample, sep = "_")),
      obs_col = c("row", "col")
    )
  }
}

sample2patient <- colData(spe)[, c("sample_id", "subject")] %>%
  as.data.frame() %>%
  dplyr::distinct() %>%
  tibble::deframe()

img_links <- c(
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151507_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151508_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151509_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151510_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151669_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151670_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151671_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151672_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151673_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151674_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151675_full_image.tif",
  "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151676_full_image.tif"
)

img_links <- tibble::as_tibble(list("link" = img_links)) %>%
  dplyr::mutate(
    sample = stringr::str_extract(link, "([^/]+)_full_image.tif$", group = 1),
    patient = sample2patient[sample],
    filename = "H_E.tiff"
  )

options(timeout = 60 * 60)
purrr::pwalk(img_links, function(link, sample, patient, filename) {
  download.file(
    link,
    file.path(out_dir, paste(patient, sample, sep = "_"), filename),
    "wget",
    quiet = TRUE
  )
})

purrr::pwalk(img_links, function(link, sample, patient, filename) {
  json <- file(file.path(out_dir, paste(patient, sample, sep = "_"), "H_E.json"))
  writeLines(c('{"scale": 1}'), json)
  close(json)
})

colData(spe) %>%
  as.data.frame() %>%
  dplyr::select(patient = subject, sample = sample_id, position, replicate, label) %>%
  dplyr::filter(!is.na(label)) %>%
  dplyr::distinct() %>%
  dplyr::count(patient, sample, position, replicate) %>%
  dplyr::rename(n_clusters = n) %>%
  dplyr::mutate(directory = paste(patient, sample, sep = "_")) %>%
  `row.names<-`(NULL) %>%
  write.table(file.path(out_dir, "samples.tsv"), sep = "\t", col.names = NA, quote = FALSE)

json <- file(file.path(out_dir, "experiment.json"))
writeLines(c('{"technology": "Visium"}'), json)
close(json)
