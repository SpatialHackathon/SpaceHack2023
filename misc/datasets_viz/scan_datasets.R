# 
# Author: Mark Robinson 
# Purpose: basically loop through all the datasets we have 
# and make some plots to give a bit of a status of what
# we are dealing with
#
# assumes packages are install (can be run with `R CMD BATCH`)
#

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
  library(Matrix)
  library(tibble)
  library(data.table)
})

source("utils.R")

# get datasets
dirs <- c("~/scratch/SpaceHack2/data",
          "~/data/results")

# make table
files_dirs <- lapply(dirs, function(u) {
  lf <- list.files(u, "counts.mtx", recursive = TRUE, full.names = TRUE)
  lf <- lf[!grepl("log1p/counts.mtx", lf, fixed = TRUE)]
  ss <- strsplit(lf, "/")
  last2 <- lapply(strsplit(lf, "/"), function(u) u %>% tail(3) %>% head(2))
  datadir <- lapply(strsplit(lf, "/"), function(u) u %>% head(length(u)-3))
  data.frame(dataset = sapply(last2, .subset, 1),
             sample = sapply(last2, .subset, 2),
             dir =  sapply(datadir, paste, collapse="/"),
             path = dirname(lf))
}) %>% 
  bind_rows %>%
  mutate(dataset_orig = dataset,
         sample = tolower(sample),
         dataset = tolower(dataset)) %>%
  arrange(sample, dataset, dir) #%>% filter(grepl("xenium",dataset_orig))


# quick tabular summary
files_dirs %>% group_by(dataset_orig, dir) %>% 
  tally() %>% arrange(dataset_orig, dir) %>% 
  #filter(grepl("xenium",dataset_orig)) %>%
  #tally() %>% arrange(desc(n)) %>% 
  print(n = Inf)

# datasets <- split(files_dirs, 
#                   paste0(files_dirs$dir, ":::", files_dirs$dataset))
datasets <- split(files_dirs, files_dirs$sample)
#datasets <- split(files_dirs, files_dirs$sample)

# colour palette with up to 15 groups
# attribution to:
# https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
 "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
 "#920000","#924900","#db6d00","#24ff24","#ffff6d")

plotlist <- list()   
count <- 0

# for(i in 1:5) {
for(i in 1:length(datasets)) {
  this_df <- datasets[[i]]
  this_df$status <- "unknown"
  annof <- TRUE
  for(j in 1:nrow(this_df)) {
    ddir <- this_df$path[j]
    cat("---------\nWorking on", ddir, "\n")

    if(!file.exists(file.path(ddir, "labels.tsv")))  annof <- FALSE

    # ----------------------
    # some extra testing
    # ----------------------
    tc <- tryCatch(readMM(file.path(ddir, "counts.mtx")), 
                   error = function(e) e)
    if( "error" %in% class(tc) ) {
        print(tc)
        this_df$status[j] <- "error-counts.mtx"
        datasets[[i]] <- this_df
        next
    }

    #tc <- tryCatch(read.delim(file.path(ddir, "observations.tsv"), 
    #                          stringsAsFactors=FALSE, row.names = 1), 
    #               error = function(e) e)
    #if( "error" %in% class(tc) ) {
    #    print(tc)
    #    this_df$status[j] <- "error-observations.tsv"
    #    datasets[[i]] <- this_df
    #    next
    #}
                   
    # read data, annotations
    # ----------------------
    spe <- get_SpatialExperiment(feature_file = file.path(ddir, "features.tsv"),
                                 matrix_file = file.path(ddir, "counts.mtx"),
                                 coord_file = file.path(ddir, "coordinates.tsv"),
                                 observation_file = file.path(ddir, "observations.tsv"))
    

    # some extra faffing
    if( is.null(rownames(colData(spe))) )
        rownames(colData(spe)) <- 1:nrow(colData(spe))
    colnames(spatialCoords(spe)) <- colnames(spatialCoords(spe)) %>% tolower

    # ----------------------
    # organize into data.frame
    # ----------------------
    df <- data.frame(rowname = rownames(colData(spe)), 
                     spatialCoords(spe),
                     total_count = colSums(counts(spe)))

    # ----------------------
    # massage labels into max 14+1 categories
    # ----------------------
    if(annof) {
        #anno <- read.delim(file.path(ddir, "labels.tsv"),
        #                   sep="\t", row.names=1)

        anno <- fread(file.path(ddir, "labels.tsv")) %>% as_tibble %>% 
                  filter(!is.na(V1)) %>%
                  column_to_rownames("V1") %>% as.data.frame
        colnames(anno)[1] <- "label"
        
        anno$rowname <- rownames(anno)
        df <- df %>% left_join(anno, by = "rowname")
        tb <- table(df$label) %>% sort(decreasing=TRUE)
        lv <- head(names(tb), min(14,length(tb)))
        if(length(tb) > 14) {
            lv <- c(lv, paste0("other_", length(tb)-14))
            df$label[df$label %in% setdiff(names(tb),lv)] <- last(lv)
        }
        df$label <- factor(df$label, lv)
        print(table(df$label, useNA="ifany"))

        this_pal <- head(pal, length(lv)) %>% setNames(lv)
        this_df$status[j] <- "good-wlabels"

    } else {
        cat("No labels.\n")
        df$label = "none"
        this_pal <- tail(pal, 1) %>% setNames("none")
        this_df$status[j] <- "good-nolabels"
    }
 
    # ----------------------
    # make a plot
    # ----------------------
    p <- ggplot(df, aes(x = x, y = y, colour = label)) + 
      geom_point() +
      theme(legend.position="bottom") +
      scale_colour_manual(values = this_pal) +
      ggtitle(paste0(this_df$dataset[j], " // ", this_df$sample[j])) +
      theme_classic()
    q <- ggplot(df, aes(x = x, y = y, colour = total_count)) + 
      geom_point() +
      theme(legend.position="bottom") + 
      ggtitle(ddir) +
      scale_colour_gradient(low = "lightgray", high = "navyblue") +
      theme_classic()

    count <- count+1
    plotlist[[count]] <- plot_grid(p, q, align = "h")

    datasets[[i]] <- this_df

  }
}

# ----------------------
# make summary table
# ----------------------
summary_tab <- bind_rows(datasets)
saveRDS(summary_tab, "summary_tab.rds")
                   
write.csv(summary_tab,
          file=paste0("spacehack-annotations_",gsub(" ","-",date()),".csv"),
          row.names = FALSE, quote = FALSE)


# ----------------------
# put all plots together
# ----------------------
pdf(paste0("spacehack-annotations_",gsub(" ","-",date()),".pdf"), 
    onefile = TRUE, width = 14, height = 6)
for (i in 1:length(plotlist))
  show(plotlist[[i]])
dev.off()



# ----------------------
# capture versions
# ----------------------
sessionInfo()


