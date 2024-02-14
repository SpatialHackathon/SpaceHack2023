#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Florian Heyl (@heylf); created code

# H_E.json and H_E.tiff not public. Request for access is still unanswered.

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(SpatialExperiment))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(WeberDivechaLCdata))
suppressPackageStartupMessages(library(Matrix))

option_list <- list(
  make_option(
    c("-o", "--out_dir"),
    type = "character", default = NULL,
    help = "Output directory to write files to."
  )
)

description <- "Load data (Visium) for Locus_coeruleus from Lukas M. Weber at al. (2022); 
The gene expression landscape of the human locus coeruleus revealed 
by single-nucleus and spatially-resolved transcriptomics"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

out_dir <- opt$out_dir

args = commandArgs(trailingOnly=TRUE)

# Load data
spe <- WeberDivechaLCdata_Visium()

coords <- as.data.frame(spe@int_colData$spatialCoords)
colnames(coords) <- c('x', 'y')
coords_rownames <- rownames(spe@int_colData$spatialCoords)

counts <- spe@assays@data$counts

counts_func <- function(x){
  fields <- unlist(strsplit(x, "_"))
  if ( length(fields) == 4 ){
    return(paste(fields[1:3], collapse = "_"))
  } else {
    return(paste(fields[1:4], collapse = "_"))
  }
}
counts_lc <- unlist(lapply(colnames(counts), counts_func))

LC_samples <- unique(spe@colData$sample_id)

for ( dir in LC_samples ){
    dir <- paste0(out_dir, "/", dir)
    if ( dir.exists(dir) == FALSE ){
        dir.create(dir, recursive=TRUE)
    }
}

print("Create output ...")

patient_list <- c()
sample_list <- c()
directory_list <- c()

# Write coordinates.tsv, observations.tsv, features.tsv, counts.mtx and labels.tsv
for (lc in LC_samples){
    
    print(lc)
    dir <- paste0(out_dir, "/", lc)

    # Write coordinates.tsv
    coords_subset <- coords[which(spe@colData$sample_id == lc),]
    rownames(coords_subset) <- coords_rownames[which(spe@colData$sample_id == lc)]
    write.table(coords_subset, file = paste0(dir, "/coordinates.tsv"), col.names = NA,
                sep = "\t", quote = FALSE, row.names = TRUE)
  
    # Count matrix has rows = genes/features and cols = cells/observations
    counts_subset <- counts[,which(counts_lc == lc)]
    
    # Transpose to have rows = cells/observations
    counts_subset <- t(counts_subset)
    writeMM(counts_subset, file = paste0(dir, "/counts.mtx"))
  
    observations_subset <- spe@colData[which(counts_lc == lc),]
    rownames(observations_subset) <- lapply(rownames(observations_subset), function(x){tail(unlist(strsplit(x,"_")),1)})
    write.table(observations_subset, file = paste0(dir, "/observations.tsv"), col.names = NA, sep = "\t", quote = FALSE)

    labels <- spe@colData$annot_region[which(spe@colData$sample_id == lc)]
    labels[which(labels == TRUE)] <- "LC"
    labels[which(labels == FALSE)] <- "non_LC"

    labels_df <- data.frame(label = labels)
    rownames(labels_df) <- rownames(observations_subset)
    write.table(labels_df, file = paste0(dir, "/labels.tsv"), col.names = NA, sep = "\t", quote = FALSE)

    # Fill metadata
    patient_list <- c(patient_list, as.character(unique(spe@colData[which(spe@colData$sample_id == lc),]$donor_id)))
    sample_list <- c(sample_list, lc)
    directory_list <- c(directory_list, dir)

    # Write features.tsv
    features <- as.data.frame(spe@rowRanges)
    rownames(features) <- spe@rowRanges$gene_id
    write.table(features, file = paste0(dir,"/features.tsv"), col.names = NA, sep = "\t", quote = FALSE)
}

## Metadata files
samples_df <- data.frame(
  patient = patient_list, 
  sample = sample_list, 
  position = rep(NA, length(patient_list)), # Not sure what position means
  replicate = rep(NA, length(patient_list)), # If they have replicated then it is really badly named
  directory = directory_list, 
  n_clusters = rep(2, length(patient_list)),
  stringsAsFactors = FALSE
)
row.names(samples_df) <- NULL
write.table(samples_df, file = file.path(out_dir, "samples.tsv"), sep = "\t", col.names = NA, quote = FALSE)

technology = "Visium"
json <- file(file.path(out_dir, "experiment.json"))
writeLines(c(paste0('{"technology": "', technology, '"}')), json)
close(json)

print("...finished")
