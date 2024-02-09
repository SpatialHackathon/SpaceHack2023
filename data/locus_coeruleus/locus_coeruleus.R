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

counts <- spe@assays@data$counts
counts_lc <- unlist(lapply(colnames(counts), function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")}))
observations_lc <- unlist(lapply(colnames(counts), function(x){paste(unlist(strsplit(x, "_"))[4], collapse="_")}))

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
    # Please check this again beause I had to flip the coordinates to get the pictures as seen in:
    # https://libd.shinyapps.io/locus-c_Visium/. See ggplot below.
    coords_subset <- coords[which(spe@colData$sample_id == lc),]
    write.table(coords_subset, file = paste0(dir, "/coordinates.tsv"), col.names = NA,
                sep = "\t", quote = FALSE, row.names = TRUE)
    
    # Create the scatter plot
    # plot_df <- coords[which(spe@colData$sample_id == lc),]
    # plot_df$label <- labels
    # plot_df
    # ggplot(plot_df, aes(x = x, y = y, color=label)) +
    #   geom_point() +
    #   coord_fixed() +  # Keep aspect ratio equal
    #   scale_y_reverse()  # Flip the y-axis

    # Count matrix has rows = genes/features and cols = cells/observations
    # Has the header %%MatrixMarket matrix coordinate integer general 
    counts_subset <- counts[,which(counts_lc == lc)]
    if (lc == "Br6522_LC_1_round1" || lc == "Br6522_LC_2_round1"){
        split_lc <- unlist(strsplit(as.character(LC_samples[1]), "_"))[1:3]
        correct_lc <- paste(split_lc, collapse = "_")
        counts_subset <- counts[,which(counts_lc == correct_lc)]
    }
    t(counts_subset)
    writeMM(counts_subset, file = paste0(dir, "/counts.mtx"))
  
    # Write observations.tsv
    observations_subset <- observations_lc[which(counts_lc == lc)]
    write.table(observations_subset, file = paste0(dir, "/observations.tsv"), col.names = NA, sep = "\t", quote = FALSE)

    labels <- spe@colData$annot_region[which(spe@colData$sample_id == lc)]
    labels[which(labels == TRUE)] <- "LC"
    labels[which(labels == FALSE)] <- "non_LC"
    write.table(labels, file = paste0(dir, "/labels.tsv"), col.names = NA, sep = "\t", quote = FALSE)

    # Fill metadata
    patient_list <- c(patient_list, as.character(unique(spe@colData[which(spe@colData$sample_id == lc),]$donor_id)))
    sample_list <- c(sample_list, lc)
    directory_list <- c(directory_list, dir)

    # Write features.tsv
    features <- rownames(spe@assays@data$counts) 
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
