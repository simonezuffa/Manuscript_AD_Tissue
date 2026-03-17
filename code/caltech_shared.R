setwd("~/Desktop/Manuscript_AD_Tissue")

library(tidyverse)
library(ggpubr)
library(patchwork)
library(rstatix)
library(UpSetR)
library(readxl)

# Read VIP tables
brain_3xtg <- read_excel("data/mice_caltech/supplementary_tables.xlsx", sheet = 1)
brain_5xfad <- read_excel("data/mice_caltech/supplementary_tables.xlsx", sheet = 7)

serum_3xtg <- read_excel("data/mice_caltech/supplementary_tables.xlsx", sheet = 2)
serum_5xfad <- read_excel("data/mice_caltech/supplementary_tables.xlsx", sheet = 8)

liver_3xtg <- read_excel("data/mice_caltech/supplementary_tables.xlsx", sheet = 3)
liver_5xfad <- read_excel("data/mice_caltech/supplementary_tables.xlsx", sheet = 9)

cecum_3xtg <- read_excel("data/mice_caltech/supplementary_tables.xlsx", sheet = 4)
cecum_5xfad <- read_excel("data/mice_caltech/supplementary_tables.xlsx", sheet = 10)

###############
# Upset plots #
###############

# Brain
list_brain <- list(
  `3xTg Up` = (brain_3xtg %>% dplyr::filter(SPF_genotype == "Mut"))$Feature_ID,
  `3xTg Down` = (brain_3xtg %>% dplyr::filter(SPF_genotype == "WT"))$Feature_ID,
  `5xFAD Up` = (brain_5xfad %>% dplyr::filter(SPF_genotype == "Mut"))$Features_ID,
  `5xFAD Down` =  (brain_5xfad %>% dplyr::filter(SPF_genotype == "WT"))$Features_ID)

brain_upset <- UpSetR::upset(fromList(list_brain), nintersects = NA,
                             point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                             sets = c("3xTg Up", "3xTg Down", "5xFAD Up", "5xFAD Down"),
                             queries = list(list(query = intersects, params = list("3xTg Up", 
                                                                                   "5xFAD Up"), color = "#287DAB", active = T), 
                                            list(query = intersects, params = list("3xTg Down", 
                                                                                   "5xFAD Down"),color = "#E5BF86", active = T)))

pdf(file="upset_brain.pdf", width = 3, height = 3)
brain_upset
dev.off()

# Serum
list_serum <- list(
  `3xTg Up` = (serum_3xtg %>% dplyr::filter(SPF_genotype == "Mut"))$Feature_ID,
  `3xTg Down` = (serum_3xtg %>% dplyr::filter(SPF_genotype == "WT"))$Feature_ID,
  `5xFAD Up` = (serum_5xfad %>% dplyr::filter(SPF_genotype == "Mut"))$Feature_ID,
  `5xFAD Down` =  (serum_5xfad %>% dplyr::filter(SPF_genotype == "WT"))$Feature_ID)

serum_upset <- UpSetR::upset(fromList(list_serum), nintersects = NA,
                             point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                             sets = c("3xTg Up", "3xTg Down", "5xFAD Up", "5xFAD Down"),
                             queries = list(list(query = intersects, params = list("3xTg Up", 
                                                                                   "5xFAD Up"), color = "#287DAB", active = T), 
                                            list(query = intersects, params = list("3xTg Down", 
                                                                                   "5xFAD Down"),color = "#E5BF86", active = T)))

pdf(file="upset_serum.pdf", width = 3, height = 3)
serum_upset
dev.off()

# Liver
list_liver <- list(
  `3xTg Up` = (liver_3xtg %>% dplyr::filter(SPF_genotype == "Mut"))$Feature_ID,
  `3xTg Down` = (liver_3xtg %>% dplyr::filter(SPF_genotype == "WT"))$Feature_ID,
  `5xFAD Up` = (liver_5xfad %>% dplyr::filter(SPF_genotype == "Mut"))$Feature_ID,
  `5xFAD Down` =  (liver_5xfad %>% dplyr::filter(SPF_genotype == "WT"))$Feature_ID)

liver_upset <- UpSetR::upset(fromList(list_liver), nintersects = NA,
                             point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                             sets = c("3xTg Up", "3xTg Down", "5xFAD Up", "5xFAD Down"),
                             queries = list(list(query = intersects, params = list("3xTg Up", 
                                                                                   "5xFAD Up"), color = "#287DAB", active = T), 
                                            list(query = intersects, params = list("3xTg Down", 
                                                                                   "5xFAD Down"),color = "#E5BF86", active = T)))

pdf(file="upset_liver.pdf", width = 3, height = 3)
liver_upset
dev.off()

# Cecum
list_cecum <- list(
  `3xTg Up` = (cecum_3xtg %>% dplyr::filter(SPF_genotype == "Mut"))$Feature_ID,
  `3xTg Down` = (cecum_3xtg %>% dplyr::filter(SPF_genotype == "WT"))$Feature_ID,
  `5xFAD Up` = (cecum_5xfad %>% dplyr::filter(SPF_genotype == "Mut"))$Feature_ID,
  `5xFAD Down` =  (cecum_5xfad %>% dplyr::filter(SPF_genotype == "WT"))$Feature_ID)

cecum_upset <- UpSetR::upset(fromList(list_cecum), nintersects = NA,
                             point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                             sets = c("3xTg Up", "3xTg Down", "5xFAD Up", "5xFAD Down"),
                             queries = list(list(query = intersects, params = list("3xTg Up", 
                                                                                   "5xFAD Up"), color = "#287DAB", active = T), 
                                            list(query = intersects, params = list("3xTg Down", 
                                                                                   "5xFAD Down"),color = "#E5BF86", active = T)))

pdf(file="upset_cecum.pdf", width = 3, height = 3)
cecum_upset
dev.off()

