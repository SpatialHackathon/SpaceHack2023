
library(readr)

summary_tab <- readRDS("summary_tab.rds")
summary_tab$status <- factor(summary_tab$status, 
                             levels=c("good-wlabels",
                                      "good-nolabels",
                                      "error-counts.mtx"))

fn <- "SpaceHack2.0 - Spatial Transcriptomics Domain Identification - Datasets.csv"

#> colnames(sh_tab)
# [1] "Technology"              "Organism"                "Tissue"                 
# [4] "state"                   "Annotation"              "Tissue\nsamples"        
# [7] "Sections \nper sample"   "Total samples"           "Histology"              
#[10] "Usable \n(GT available)" "EULA or \nLicense.."     "Citation/Reference"     
#[13] "Link"                    "folder"                  "Notes" 

sh_tab <- read_csv(fn) %>%
    select(Technology, Organism, Tissue, state, Annotation, folder)

x <- summary_tab %>% 
    group_by(dataset_orig, dir) %>% 
    summarize(status=paste0(table(status), collapse=";")) %>%
    filter(row_number()==1) 

z <- x %>%
    left_join(sh_tab, by = c("dataset_orig" = "folder")) %>%
    data.frame

z %>% filter(is.na(Organism)) %>% 
    select(dataset_orig, dir, status, state, Technology) %>%
    arrange(desc(status))

z %>% filter(!is.na(Organism)) %>% 
     select(dataset_orig, dir, status, state, Technology) %>%
     mutate(state=substr(state, 1, 10))

sh_tab %>% left_join(x, by = c("folder" = "dataset_orig")) %>% 
    filter(is.na(folder)) %>% arrange(desc(state)) %>%
    select(-dir, -status, -folder)

keep <- z %>% filter(!is.na(Organism), grepl("^[1-9]", status)) %>% 
    select(dataset_orig, dir, status, state, Technology) %>% pull(dataset_orig)

saveRDS(keep, "keep.rds")



