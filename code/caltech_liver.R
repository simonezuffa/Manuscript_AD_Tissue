setwd("~/OneDrive - University of California, San Diego Health/Projects/Caltech/Manuscript_AD_Tissue")

library(tidyverse)
library(mixOmics)
library(ggpubr)
library(vegan)
library(caret)
library(limma)
library(patchwork)
library(rstatix)
library(Spectra)
library(MsBackendMgf)
library(UpSetR)
library(effsize)

# Read in data
feature_table <- read_csv("data/mice_caltech/liver/gnps_quant_liver.csv")
metadata <- read_csv("data/mice_caltech/liver/meta_liver.csv")
sample_order <- read.csv("data/mice_caltech/liver/sequence_liver.csv")
annotations <- read.delim("data/mice_caltech/liver/fbmn/nf_output/library/merged_results_with_gnps.tsv") %>%
  dplyr::filter(!str_detect(pattern = "REFRAME", LibraryName)) # remove REFRAME library
annotations$X.Scan. <- as.character(annotations$X.Scan.)
canopus <- read_tsv("data/mice_caltech/liver/canopus_liver.tsv") %>% dplyr::select(1,3,5,7,9)
canopus$id <- gsub("^.*sirius_", "", canopus$id)

sample_order$Plate <- gsub("^(.*?):.*", "\\1", sample_order$Position)

info_feature <- feature_table %>% dplyr::select(1:3,7)
colnames(info_feature) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature$Feature <- as.character(info_feature$Feature)

info_feature_complete <- info_feature %>% 
  left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  dplyr::select(1:5,18,24) %>% 
  left_join(canopus %>% distinct(id, .keep_all = TRUE), by = c("Feature" = "id"))

#write_csv(x = info_feature_complete, file = "liver_annotations.csv")

# Extra bile acids libraries and rev cosine
ba_ipsita <- read_csv("data/mice_caltech/liver/merged_Bile_acid_classic_networking.csv") %>% 
  dplyr::filter(SpectrumID %in% annotations$SpectrumID)
ba_rev <- read_tsv("data/mice_caltech/liver/candidate_BA_with_best_annotations.tsv") %>% 
  dplyr::filter(query_id %in% annotations$SpectrumID)


# Data table
data <- feature_table %>%
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE)

data$SampleID <- gsub(".mzML Peak area", "", data$SampleID)

# Metadata
metadata_metabolomics <- data_frame(SampleID = data$SampleID) %>% 
  left_join(metadata, by = c("SampleID" = "Alias")) %>% 
  left_join(sample_order, by = c("SampleID" = "File.Name")) %>% 
  dplyr::select(-c("sample_type", "intervention")) %>%
  dplyr::mutate(Colonized = case_when(str_detect(pattern = "GF", mouse_strain) ~ "GF",
                                      str_detect(pattern = "SPF", mouse_strain) ~ "SPF",
                                      str_detect(pattern = "Blank_", SampleID) ~ "Blank",
                                      str_detect(pattern = "Pool_qc_", SampleID) ~ "Pool",
                                      str_detect(pattern = "sixmix_", SampleID) ~ "SixMix",
                                      str_detect(pattern = "srm_", SampleID) ~ "SRM",
                                      TRUE ~ "blank/qc")) %>%
  dplyr::mutate(Study = case_when(str_detect(pattern = "GF HET|SPF HET|GF WT|SPF WT", mouse_strain) ~ "5xFAD",
                                  str_detect(pattern = "GF 3XTG|SPF 3XTG|GF B6|SPF B6", mouse_strain) ~ "3xTG",
                                  str_detect(pattern = "Blank_", SampleID) ~ "Blank",
                                  str_detect(pattern = "Pool_qc_", SampleID) ~ "Pool",
                                  str_detect(pattern = "sixmix_", SampleID) ~ "SixMix",
                                  str_detect(pattern = "srm_", SampleID) ~ "SRM",
                                  TRUE ~ "blank/qc")) %>%
  dplyr::mutate(Strain = case_when(str_detect(pattern = "WT|B6", mouse_strain) ~ "WT",
                                   str_detect(pattern = "3XTG|HET", mouse_strain) ~ "Mut",
                                   str_detect(pattern = "Blank_", SampleID) ~ "Blank",
                                   str_detect(pattern = "Pool_qc_", SampleID) ~ "Pool",
                                   str_detect(pattern = "sixmix_", SampleID) ~ "SixMix",
                                   str_detect(pattern = "srm_", SampleID) ~ "SRM",
                                   TRUE ~ "blank/qc"))


# Investigate total peak area
data_TIC <- data.frame(TIC = rowSums(data %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("File.Name") %>% left_join(sample_order)

data_TIC %>% ggscatter("Order", "TIC", add = "reg.line") +
  stat_cor()

# Problem in the middle of the run --> samples have been re-run

# Sample rerun
sample_rerun <- data %>% dplyr::filter(str_detect(pattern = "_2024", SampleID)) %>%
  dplyr::select(SampleID)

sample_rerun$SampleID1 <- sub("_[^_]*$", "", sample_rerun$SampleID)

data_filter <- data %>% dplyr::filter(!(SampleID %in% sample_rerun$SampleID1))
data_filter$SampleID <- sub("_2024.*", "", data_filter$SampleID)

# Investigate total peak area
data_TIC <- data.frame(TIC = rowSums(data_filter %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("File.Name") %>% left_join(sample_order)

data_TIC %>% dplyr::filter(str_detect(pattern = "3xtg|5xfad", File.Name)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 2e9) +
  stat_cor()


# Check sample type
sample_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "3xtg|5xfad", File.Name)) %>% summarise(mean(TIC))
pool_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Pool", File.Name)) %>% summarise(mean(TIC))
six_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "sixmix", File.Name)) %>% summarise(mean(TIC))
blank_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", File.Name)) %>% summarise(mean(TIC))
SRM_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "srm", File.Name)) %>% summarise(mean(TIC))

ratio_tic_pb <- pool_tic/sample_tic
ratio_tic_sb <- sample_tic/blank_tic
ratio_tic_6p <- pool_tic/six_tic

# Check TIC overtime for QCpool, QCmix, Blank, SRM
data_TIC %>% dplyr::filter(str_detect(pattern = "Pool", File.Name)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 2e9) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "sixmix", File.Name)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 1e9) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", File.Name)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 7e8) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "srm", File.Name)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 1e9) +
  stat_cor()


# Check internal standard and remove sample where there was a shift
fbmn_IS <- annotations %>% dplyr::filter(str_detect(Compound_Name, regex("sulf", ignore_case = TRUE))) %>% 
  distinct(X.Scan., .keep_all = TRUE) %>% dplyr::filter(Organism != "BILELIB19")

# Extract IS features for the table
table_IS <- data_filter %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% 
  dplyr::filter(ID %in% fbmn_IS$X.Scan.) %>% column_to_rownames("ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SampleID") %>% dplyr::filter(!(str_detect(SampleID, "srm|sixmix|Blank|blank"))) %>%
  dplyr::select(SampleID, `8436`) %>% left_join(metadata_metabolomics)

colnames(table_IS)[2] <- "Sulfadimethoxine"

table_IS %>% ggscatter(x = "Order", y = "Sulfadimethoxine", add = c("reg.line")) +
  ylim(0, 9e6) + stat_cor()

table_IS %>% ggbarplot(x = "Order", y = "Sulfadimethoxine", xlab = "Run Order", 
                       ylab = "Peak Area Sulfadimethoxine", title = "Internal Standard Acquisition") +
  geom_hline(yintercept = mean(table_IS$Sulfadimethoxine, na.rm = TRUE), linetype = "dashed", color = "blue")

cv_is <- sd(table_IS$Sulfadimethoxine)/mean(table_IS$Sulfadimethoxine)


# Check features per sample type
data_blank <- data_filter %>% dplyr::filter(str_detect(pattern = "Blank", SampleID))
data_pool <- data_filter %>% dplyr::filter(str_detect(pattern = "Pool", SampleID))
data_sixmix <- data_filter %>% dplyr::filter(str_detect(pattern = "sixmix", SampleID))
data_srm <- data_filter %>% dplyr::filter(str_detect(pattern = "srm", SampleID))

# Blank
blanks_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                  Mean_blank = data_blank %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_blank =  data_blank %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_blank, SD_blank, CV_blank) %>% 
  dplyr::filter(Mean_blank > 0) %>% arrange(desc(Mean_blank))

# Six mix
sixmix_feature_info <- data.frame(Feature = colnames(data_sixmix)[-1],
                                  Mean_sixmix = data_sixmix %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_sixmix = data_sixmix %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sixmix = SD_sixmix/Mean_sixmix) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_sixmix, SD_sixmix, CV_sixmix) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% arrange(desc(Mean_sixmix))
# CVs of the six standards are ranging between 19% to 22% -> ok

# Pool
pools_feature_info <- data.frame(Feature = colnames(data_pool)[-1],
                                 Mean_pool = data_pool %>% column_to_rownames("SampleID") %>% colMeans(), 
                                 SD_pool =  data_pool %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_pool = SD_pool/Mean_pool) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_pool, SD_pool, CV_pool) %>% 
  dplyr::filter(Mean_pool > 0) %>% arrange(desc(Mean_pool))


# Features to be removed Pools/Blank < 5
feature_to_remove <- blanks_feature_info %>% left_join(pools_feature_info) %>%
  dplyr::filter(Mean_blank > 0) %>% 
  dplyr::mutate(Pool_Blank = Mean_pool/Mean_blank) %>% 
  dplyr::filter(Pool_Blank < 5 | is.na(Pool_Blank))

# Data with blank removal
data_clean <- data_filter %>% dplyr::select(-c(feature_to_remove$Feature)) %>%
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")

# Features to be removed Pool/QCmix < 5
feature_to_remove_mix <- sixmix_feature_info %>% left_join(pools_feature_info) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% 
  dplyr::mutate(Pool_Mix = Mean_pool/Mean_sixmix) %>% 
  dplyr::filter(Pool_Mix < 5 | is.na(Pool_Mix)) %>% dplyr::filter(!(Feature %in% feature_to_remove$Feature))

# Data with QCmix removal
data_clean2 <- data_clean %>% dplyr::select(-c(feature_to_remove_mix$Feature))

# Remove feature before 0.2 minutes and after 8 minutes
feature_to_rt <- info_feature_complete %>% dplyr::filter(RT < 0.2 | RT > 8) %>%
  dplyr::filter(!(Feature %in% feature_to_remove$Feature)) %>%
  dplyr::filter(!(Feature %in% feature_to_remove_mix$Feature)) %>%
  dplyr::filter(!(Feature %in% c(16911, 16896, 16258)))

# Final cleaned table
data_clean3 <- data_clean2 %>% dplyr::select(-c(feature_to_rt$Feature)) %>%
  dplyr::filter(!(SampleID %in% c("L3xtg_156", "L3xtg_147", "L3xtg_127"))) #outliers

# PCA raw data
PCA_raw <- mixOmics::pca(data_clean3 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

i <- "Study"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Liver", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Keep only samples
data_sample <- data_clean3 %>% 
  dplyr::filter(!(str_detect(pattern = "srm|sixmix|Pool|Blank|blank", SampleID)))

data_sample <- data_clean3 %>% 
  dplyr::filter(!(str_detect(pattern = "srm|sixmix|Pool|Blank|blank", SampleID))) %>%
  dplyr::filter(!(SampleID %in% c("L5xfad_1", "L5xfad_10", "L5xfad_11", "L5xfad_12", 
                                  "L5xfad_13", "L5xfad_14", "L5xfad_15", "L5xfad_16", 
                                  "L5xfad_17", "L5xfad_18", "L5xfad_19", "L5xfad_2",
                                  "L5xfad_20", "L5xfad_21", "L5xfad_22", "L5xfad_23",
                                  "L5xfad_24", "L5xfad_25", "L5xfad_26", "L5xfad_27", 
                                  "L5xfad_28", "L5xfad_29", "L5xfad_3", "L5xfad_30",
                                  "L5xfad_4", "L5xfad_5", "L5xfad_6", "L5xfad_7",
                                  "L5xfad_8", "L5xfad_9"))) 
# this is to generate the figure for manuscript removing corrupted samples in box 9
# from the 5xFAD study (confirmed also via biocrates)

# PCA raw data
PCA_raw <- mixOmics::pca(data_sample %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

i <- "Study"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Liver ", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# RCLR transformation
data_sample_clr <- decostand(data_sample %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

i <- "Study"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = "Liver", palette = c("#009392", "#F0746E"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = PCA_plot, filename = "liver_pca_metabolome.svg", device = "svg", dpi = "retina", height = 2, width = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_sample_clr, method = "euclidean")
disper_study <- betadisper(dist_metabolites, PCA_whole_scores$Study)
anova(disper_study)
permanova <- adonis2(dist_metabolites ~ Study, PCA_whole_scores, na.action = na.omit)

#write_csv(x = data_sample, file = "liver_metabolome.csv")
#write_csv(x = data_sample_clr, file = "liver_metabolome_rclr.csv")
#write_csv(x = metadata_metabolomics %>% dplyr::filter(SampleID %in% data_sample$SampleID), file = "liver_metadata.csv")

# Work separately on the two different studies cause there is a strong clustering 

# Separate animal by study, colonization and strain
sample_3xtg <- metadata_metabolomics %>% dplyr::filter(Study == "3xTG")
sample_5xfad <- metadata_metabolomics %>% dplyr::filter(Study == "5xFAD")

sample_gf <- metadata_metabolomics %>% dplyr::filter(Colonized == "GF")
sample_spf <- metadata_metabolomics %>% dplyr::filter(Colonized == "SPF")

sample_wt <- metadata_metabolomics %>% dplyr::filter(Strain == "WT")
sample_mut <- metadata_metabolomics %>% dplyr::filter(Strain == "Mut")

sample_male <- metadata_metabolomics %>% dplyr::filter(sex == "male")
sample_female <- metadata_metabolomics %>% dplyr::filter(sex == "female")


########
# 3xTG #
########
data_3xtg <- data_sample %>% dplyr::filter(SampleID %in% sample_3xtg$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L3xtg_43", "L3xtg_55", "L3xtg_57", "L3xtg_59", 
                                  "L3xtg_49", "L3xtg_144", "L3xtg_117")))

# RCLR transformation
data_3xtg_clr <- decostand(data_3xtg %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_3xtg_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_3xtg_plots <- list()

for (i in c("box", "host_age", "sex", "Plate",  "Colonized", "Strain")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 3xTG", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_3xtg_plots[[i]] <- PCA_plot
  
}

PCA_3xtg_plots_final <- wrap_plots(PCA_3xtg_plots, nrow = 3)

# PERMANOVA
dist_metabolites <- vegdist(data_3xtg_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$mouse_strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ Strain * Colonized * sex + 
                       box + Plate, PCA_whole_scores, na.action = na.omit)


# PLSDA 3xTG - Strain
PLSDA_3xtg_strain <- mixOmics::plsda(data_3xtg_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                     PCA_whole_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_3xtg_strain_scores <- data.frame(PLSDA_3xtg_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_strain_plot <- PLSDA_3xtg_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 3xTG - Strain",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_strain <- plotLoadings(PLSDA_3xtg_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_strain <- perf(PLSDA_3xtg_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_strain, legend = FALSE)

VIPs_3xtg_strain <- as.data.frame(mixOmics::vip(PLSDA_3xtg_strain))
VIPs_3xtg_strain_filter <- dplyr::filter(VIPs_3xtg_strain, VIPs_3xtg_strain$comp1 > 1)
VIPs_3xtg_strain_filter$ID <- rownames(VIPs_3xtg_strain_filter)
VIPs_3xtg_strain_select <- VIPs_3xtg_strain_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_strain_Load <- VIPs_3xtg_strain_select %>% 
  left_join(Loadings_3xtg_strain, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 3xTG - Colonized
PLSDA_3xtg_colonized <- mixOmics::plsda(data_3xtg_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                        PCA_whole_scores$Colonized, ncomp = 3, scale = TRUE)
PLSDA_3xtg_colonized_scores <- data.frame(PLSDA_3xtg_colonized$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_colonized_plot <- PLSDA_3xtg_colonized_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Colonized", alpha = 0.6, title = "PLSDA - Liver 3xTG - Colonized",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_colonized$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_colonized$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_colonized_scores %>% group_by(Colonized) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Colonized), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_colonized <- plotLoadings(PLSDA_3xtg_colonized, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_colonized <- perf(PLSDA_3xtg_colonized, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_colonized, legend = FALSE)

VIPs_3xtg_colo <- as.data.frame(mixOmics::vip(PLSDA_3xtg_colonized))
VIPs_3xtg_colo_filter <- dplyr::filter(VIPs_3xtg_colo, VIPs_3xtg_colo$comp1 > 1)
VIPs_3xtg_colo_filter$ID <- rownames(VIPs_3xtg_colo_filter)
VIPs_3xtg_colo_select <- VIPs_3xtg_colo_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_colo_Load <- VIPs_3xtg_colo_select %>% 
  left_join(Loadings_3xtg_colonized, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 3xTG - Sex
PLSDA_3xtg_sex <- mixOmics::plsda(data_3xtg_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                  PCA_whole_scores$sex, ncomp = 3, scale = TRUE)
PLSDA_3xtg_sex_scores <- data.frame(PLSDA_3xtg_sex$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_sex_plot <- PLSDA_3xtg_sex_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "sex", alpha = 0.6, title = "PLSDA - Liver 3xTG - Sex",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_sex$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_sex$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_sex_scores %>% group_by(sex) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = sex), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_sex <- plotLoadings(PLSDA_3xtg_sex, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_sex <- perf(PLSDA_3xtg_sex, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_sex, legend = FALSE)

VIPs_3xtg_sex <- as.data.frame(mixOmics::vip(PLSDA_3xtg_sex))
VIPs_3xtg_sex_filter <- dplyr::filter(VIPs_3xtg_sex, VIPs_3xtg_sex$comp1 > 1)
VIPs_3xtg_sex_filter$ID <- rownames(VIPs_3xtg_sex_filter)
VIPs_3xtg_sex_select <- VIPs_3xtg_sex_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_sex_Load <- VIPs_3xtg_sex_select %>% 
  left_join(Loadings_3xtg_sex, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))


# Remove feature influenced by sex
VIPs_3xtg_strain_no_sex <- VIPs_3xtg_strain_Load %>% 
  anti_join(VIPs_3xtg_sex_Load %>% dplyr::select(ID))
VIPs_3xtg_colo_no_sex <- VIPs_3xtg_colo_Load %>% 
  anti_join(VIPs_3xtg_sex_Load %>% dplyr::select(ID))

VIPs_3xtg_colo_strain_clean <- VIPs_3xtg_strain_no_sex %>% 
  inner_join(VIPs_3xtg_colo_no_sex %>% dplyr::select(ID))

#write_csv(x = VIPs_3xtg_strain_no_sex, file = "3xtg_liver_plsda_vip_all_samples_genotype.csv")
#write_csv(x = VIPs_3xtg_colo_no_sex, file = "3xtg_liver_plsda_vip_all_samples_colonization.csv")
#write_csv(x = VIPs_3xtg_colo_strain_clean, file = "3xtg_liver_plsda_vip_all_samples_genotype_colonization_intersection.csv")


############
# Compare 3xtg gf v wt gf to 3xtg spf v wt spf to see genotype effect
############
data_3xtg_gf <- data_sample %>% dplyr::filter(SampleID %in% sample_3xtg$SampleID) %>%
  dplyr::filter(SampleID %in% sample_gf$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L3xtg_43", "L3xtg_55", "L3xtg_57", "L3xtg_59", 
                                  "L3xtg_49", "L3xtg_144", "L3xtg_117")))

data_3xtg_spf <- data_sample %>% dplyr::filter(SampleID %in% sample_3xtg$SampleID) %>%
  dplyr::filter(SampleID %in% sample_spf$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L3xtg_43", "L3xtg_55", "L3xtg_57", "L3xtg_59", 
                                  "L3xtg_49", "L3xtg_144", "L3xtg_117")))

# RCLR transformation
data_3xtg_gf_clr <- decostand(data_3xtg_gf %>% column_to_rownames("SampleID"), method = "rclr")
data_3xtg_spf_clr <- decostand(data_3xtg_spf %>% column_to_rownames("SampleID"), method = "rclr")

# PCA GF
PCA_3xtg_gf <- mixOmics::pca(data_3xtg_gf_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                             ncomp = 2, center = TRUE, scale = TRUE)
PCA_3xtg_gf_scores <- data.frame(PCA_3xtg_gf$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_3xtg_gf_plots <- list()

for (i in c("box", "host_age", "sex", "Plate", "Strain")) {
  
  PCA_plot <- PCA_3xtg_gf_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 3xTG - GF", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_3xtg_gf$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_3xtg_gf$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_3xtg_gf_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_3xtg_gf_plots[[i]] <- PCA_plot
  
}

PCA_3xtg_gf_plots_final <- wrap_plots(PCA_3xtg_gf_plots, nrow = 3)

# PCA SPF
PCA_3xtg_spf <- mixOmics::pca(data_3xtg_spf_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                              ncomp = 2, center = TRUE, scale = TRUE)
PCA_3xtg_spf_scores <- data.frame(PCA_3xtg_spf$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_3xtg_spf_plots <- list()

for (i in c("box", "host_age", "sex", "Plate", "Strain")) {
  
  PCA_plot <- PCA_3xtg_spf_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 3xTG SPF", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_3xtg_spf$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_3xtg_spf$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_3xtg_spf_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_3xtg_spf_plots[[i]] <- PCA_plot
  
}

PCA_3xtg_spf_plots_final <- wrap_plots(PCA_3xtg_spf_plots, nrow = 3)


# PLSDA 3xTG GF - Strain
PLSDA_3xtg_gf_strain <- mixOmics::plsda(data_3xtg_gf_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                        PCA_3xtg_gf_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_3xtg_gf_strain_scores <- data.frame(PLSDA_3xtg_gf_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_gf_strain_plot <- PLSDA_3xtg_gf_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 3xTG GF - Strain",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_gf_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_gf_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_gf_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_gf_strain <- plotLoadings(PLSDA_3xtg_gf_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_gf_strain <- perf(PLSDA_3xtg_gf_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_gf_strain, legend = FALSE)

VIPs_3xtg_gf_strain <- as.data.frame(mixOmics::vip(PLSDA_3xtg_gf_strain))
VIPs_3xtg_gf_strain_filter <- dplyr::filter(VIPs_3xtg_gf_strain, VIPs_3xtg_gf_strain$comp1 > 1)
VIPs_3xtg_gf_strain_filter$ID <- rownames(VIPs_3xtg_gf_strain_filter)
VIPs_3xtg_gf_strain_select <- VIPs_3xtg_gf_strain_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_gf_strain_Load <- VIPs_3xtg_gf_strain_select %>% 
  left_join(Loadings_3xtg_gf_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 3xTG SPF - Strain
PLSDA_3xtg_spf_strain <- mixOmics::plsda(data_3xtg_spf_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                         PCA_3xtg_spf_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_3xtg_spf_strain_scores <- data.frame(PLSDA_3xtg_spf_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_spf_strain_plot <- PLSDA_3xtg_spf_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 3xTG SPF - Strain",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_spf_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_spf_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_spf_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_spf_strain <- plotLoadings(PLSDA_3xtg_spf_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_spf_strain <- perf(PLSDA_3xtg_spf_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_spf_strain, legend = FALSE)

#write_csv(x = perf_plsda_3xtg_spf_strain$error.rate$overall %>% 
#            as.data.frame() %>% dplyr::mutate(Model = "3xtg_spf_strain"), file = "liver_plsda_performance_3xtg.csv")


VIPs_3xtg_spf_strain <- as.data.frame(mixOmics::vip(PLSDA_3xtg_spf_strain))
VIPs_3xtg_spf_strain_filter <- dplyr::filter(VIPs_3xtg_spf_strain, VIPs_3xtg_spf_strain$comp1 > 1)
VIPs_3xtg_spf_strain_filter$ID <- rownames(VIPs_3xtg_spf_strain_filter)
VIPs_3xtg_spf_strain_select <- VIPs_3xtg_spf_strain_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_spf_strain_Load <- VIPs_3xtg_spf_strain_select %>% 
  left_join(Loadings_3xtg_spf_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# Sample matching 
list_3xtg_strain <- list(
  `GF 3xTG` = (VIPs_3xtg_gf_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `GF WT` =  (VIPs_3xtg_gf_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID,
  `SPF 3xTG` = (VIPs_3xtg_spf_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `SPF WT` = (VIPs_3xtg_spf_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID)

strain_3xtg <- UpSetR::upset(fromList(list_3xtg_strain), nsets = 4, nintersects = 6, 
                             point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                             sets = c("GF 3xTG", "GF WT", "SPF 3xTG", "SPF WT"),
                             queries = list(list(query = intersects, params = list("SPF 3xTG", 
                                                                                   "GF 3xTG"), color = "steelblue1", active = T), 
                                            list(query = intersects, params = list("SPF WT", 
                                                                                   "GF WT"),color = "orange", active = T)))

#pdf(file="3xtg_liver_upset_genotype.pdf", width = 3, height = 2.5)
#strain_3xtg
#dev.off()


# Check common features 3xtg mutation or wt
features_3xtg_mut <- VIPs_3xtg_spf_strain_Load %>% 
  dplyr::filter(GroupContrib == "Mut") %>%
  dplyr::select(ID, comp1) %>% 
  inner_join(VIPs_3xtg_gf_strain_Load %>% 
               dplyr::filter(GroupContrib == "Mut"), by = c("ID" = "ID"))

features_3xtg_wt <- VIPs_3xtg_spf_strain_Load %>% 
  dplyr::filter(GroupContrib == "WT") %>%
  dplyr::select(ID, comp1) %>% 
  inner_join(VIPs_3xtg_gf_strain_Load %>% 
               dplyr::filter(GroupContrib == "WT"), by = c("ID" = "ID"))

#write_csv(x = features_3xtg_mut, file = "3xtg_liver_genotype_mut.csv")
#write_csv(x = features_3xtg_wt, file = "3xtg_liver_genotype_wt.csv")


############
# Compare 3xtg gf v 3xtg spf to wt gf v wt spf to see microbiome effect
############
data_3xtg_mut <- data_sample %>% dplyr::filter(SampleID %in% sample_3xtg$SampleID) %>%
  dplyr::filter(SampleID %in% sample_mut$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L3xtg_43", "L3xtg_55", "L3xtg_57", "L3xtg_59", 
                                  "L3xtg_49", "L3xtg_144", "L3xtg_117")))

data_3xtg_wt <- data_sample %>% dplyr::filter(SampleID %in% sample_3xtg$SampleID) %>%
  dplyr::filter(SampleID %in% sample_wt$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L3xtg_43", "L3xtg_55", "L3xtg_57", "L3xtg_59", 
                                  "L3xtg_49", "L3xtg_144", "L3xtg_117")))

# RCLR transformation
data_3xtg_mut_clr <- decostand(data_3xtg_mut %>% column_to_rownames("SampleID"), method = "rclr")
data_3xtg_wt_clr <- decostand(data_3xtg_wt %>% column_to_rownames("SampleID"), method = "rclr")

# PCA Mut
PCA_3xtg_mut <- mixOmics::pca(data_3xtg_mut_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                              ncomp = 2, center = TRUE, scale = TRUE)
PCA_3xtg_mut_scores <- data.frame(PCA_3xtg_mut$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_3xtg_mut_plots <- list()

for (i in c("box", "host_age", "Plate", "Colonized", "sex")) {
  
  PCA_plot <- PCA_3xtg_mut_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 3xTG - Mut", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_3xtg_mut$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_3xtg_mut$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_3xtg_mut_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_3xtg_mut_plots[[i]] <- PCA_plot
  
}

PCA_3xtg_mut_plots_final <- wrap_plots(PCA_3xtg_mut_plots, nrow = 3)

# PCA WT
PCA_3xtg_wt <- mixOmics::pca(data_3xtg_wt_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                             ncomp = 2, center = TRUE, scale = TRUE)
PCA_3xtg_wt_scores <- data.frame(PCA_3xtg_wt$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_3xtg_wt_plots <- list()

for (i in c("box", "host_age", "Plate", "Colonized", "sex")) {
  
  PCA_plot <- PCA_3xtg_wt_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 3xTG - WT", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_3xtg_wt$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_3xtg_wt$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_3xtg_wt_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_3xtg_wt_plots[[i]] <- PCA_plot
  
}

PCA_3xtg_wt_plots_final <- wrap_plots(PCA_3xtg_wt_plots, nrow = 3)


# PLSDA 3xTG Mut - Colonization
PLSDA_3xtg_mut_colo <- mixOmics::plsda(data_3xtg_mut_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                         PCA_3xtg_mut_scores$Colonized, ncomp = 3, scale = TRUE)
PLSDA_3xtg_mut_colo_scores <- data.frame(PLSDA_3xtg_mut_colo$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_mut_colo_plot <- PLSDA_3xtg_mut_colo_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Colonized", alpha = 0.6, title = "PLSDA - Liver 3xTG Mut - Colonization",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_mut_colo$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_mut_colo$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_mut_colo_scores %>% group_by(Colonized) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Colonized), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_mut_strain <- plotLoadings(PLSDA_3xtg_mut_colo, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_mut_colonization <- perf(PLSDA_3xtg_mut_colo, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_mut_colonization, legend = FALSE)

VIPs_3xtg_mut_colo <- as.data.frame(mixOmics::vip(PLSDA_3xtg_mut_colo))
VIPs_3xtg_mut_colo_filter <- dplyr::filter(VIPs_3xtg_mut_colo, VIPs_3xtg_mut_colo$comp1 > 1)
VIPs_3xtg_mut_colo_filter$ID <- rownames(VIPs_3xtg_mut_colo_filter)
VIPs_3xtg_mut_colo_select <- VIPs_3xtg_mut_colo_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_mut_colo_Load <- VIPs_3xtg_mut_colo_select %>% 
  left_join(Loadings_3xtg_mut_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 3xTG WT - Colonization
PLSDA_3xtg_wt_colo <- mixOmics::plsda(data_3xtg_wt_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                        PCA_3xtg_wt_scores$Colonized, ncomp = 2, scale = TRUE)
PLSDA_3xtg_wt_colo_scores <- data.frame(PLSDA_3xtg_wt_colo$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_wt_colo_plot <- PLSDA_3xtg_wt_colo_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Colonized", alpha = 0.6, title = "PLSDA - Liver 3xTG WT - Colonization",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_wt_colo$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_wt_colo$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_wt_colo_scores %>% group_by(Colonized) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Colonized), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_wt_strain <- plotLoadings(PLSDA_3xtg_wt_colo, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_wt_colonization <- perf(PLSDA_3xtg_wt_colo, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_wt_colonization, legend = FALSE)

VIPs_3xtg_wt_colo <- as.data.frame(mixOmics::vip(PLSDA_3xtg_wt_colo))
VIPs_3xtg_wt_colo_filter <- dplyr::filter(VIPs_3xtg_wt_colo, VIPs_3xtg_wt_colo$comp1 > 1)
VIPs_3xtg_wt_colo_filter$ID <- rownames(VIPs_3xtg_wt_colo_filter)
VIPs_3xtg_wt_colo_select <- VIPs_3xtg_wt_colo_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_wt_colo_Load <- VIPs_3xtg_wt_colo_select %>% 
  left_join(Loadings_3xtg_wt_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# Sample matching 
list_3xtg_microbiome <- list(
  `3xTG GF` = (VIPs_3xtg_mut_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `3xTG SPF` = (VIPs_3xtg_mut_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID,
  `WT GF` =  (VIPs_3xtg_wt_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `WT SPF` = (VIPs_3xtg_wt_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID)

colonization_3xtg <- UpSetR::upset(fromList(list_3xtg_microbiome), nsets = 4, nintersects = 6, 
                                   point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                                   sets = c("3xTG GF", "3xTG SPF", "WT GF", "WT SPF"),
                                   queries = list(list(query = intersects, params = list("WT SPF", 
                                                                                         "3xTG SPF"), color = "steelblue1", active = T), 
                                                  list(query = intersects, params = list("WT GF", 
                                                                                         "3xTG GF"),color = "orange", active = T)))

#pdf(file="3xtg_liver_upset_colonization.pdf", width = 3, height = 2.5)
#colonization_3xtg
#dev.off()

# Check common features spf or gf
features_3xtg_spf <- VIPs_3xtg_mut_colo_Load %>% 
  dplyr::filter(GroupContrib == "SPF") %>%
  dplyr::select(ID, comp1) %>% 
  inner_join(VIPs_3xtg_wt_colo_Load %>% 
               dplyr::filter(GroupContrib == "SPF"), by = c("ID" = "ID"))

features_3xtg_gf <- VIPs_3xtg_mut_colo_Load %>% 
  dplyr::filter(GroupContrib == "GF") %>%
  dplyr::select(ID, comp1) %>% 
  inner_join(VIPs_3xtg_wt_colo_Load %>% 
               dplyr::filter(GroupContrib == "GF"), by = c("ID" = "ID"))

#write_csv(x = features_3xtg_spf, file = "3xtg_liver_colonization_spf.csv")
#write_csv(x = features_3xtg_gf, file = "3xtg_liver_colonization_gf.csv")


############
# Compare 3xtg male v wt male to 3xtg female v wt female to see strain effect of sex
############
data_3xtg_male <- data_sample %>% dplyr::filter(SampleID %in% sample_3xtg$SampleID) %>%
  dplyr::filter(SampleID %in% sample_male$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L3xtg_43", "L3xtg_55", "L3xtg_57", "L3xtg_59", 
                                  "L3xtg_49", "L3xtg_144", "L3xtg_117")))

data_3xtg_female <- data_sample %>% dplyr::filter(SampleID %in% sample_3xtg$SampleID) %>%
  dplyr::filter(SampleID %in% sample_female$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L3xtg_43", "L3xtg_55", "L3xtg_57", "L3xtg_59", 
                                  "L3xtg_49", "L3xtg_144", "L3xtg_117")))

# RCLR transformation
data_3xtg_male_clr <- decostand(data_3xtg_male %>% column_to_rownames("SampleID"), method = "rclr")
data_3xtg_female_clr <- decostand(data_3xtg_female %>% column_to_rownames("SampleID"), method = "rclr")

# PCA Male
PCA_3xtg_male <- mixOmics::pca(data_3xtg_male_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                               ncomp = 2, center = TRUE, scale = TRUE)
PCA_3xtg_male_scores <- data.frame(PCA_3xtg_male$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_3xtg_male_plots <- list()

for (i in c("box", "host_age", "Plate", "Strain", "Colonized")) {
  
  PCA_plot <- PCA_3xtg_male_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 3xTG Male", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_3xtg_male$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_3xtg_male$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_3xtg_male_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_3xtg_male_plots[[i]] <- PCA_plot
  
}

PCA_3xtg_male_plots_final <- wrap_plots(PCA_3xtg_male_plots, nrow = 3)

# PCA Female
PCA_3xtg_female <- mixOmics::pca(data_3xtg_female_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                 ncomp = 2, center = TRUE, scale = TRUE)
PCA_3xtg_female_scores <- data.frame(PCA_3xtg_female$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_3xtg_female_plots <- list()

for (i in c("box", "host_age", "Plate", "Strain", "Colonized")) {
  
  PCA_plot <- PCA_3xtg_female_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 3xTG Female", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_3xtg_female$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_3xtg_female$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_3xtg_female_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_3xtg_female_plots[[i]] <- PCA_plot
  
}

PCA_3xtg_female_plots_final <- wrap_plots(PCA_3xtg_female_plots, nrow = 3)


# PLSDA 3xTG Male - Strain
PLSDA_3xtg_male_strain <- mixOmics::plsda(data_3xtg_male_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                          PCA_3xtg_male_scores$Strain, ncomp = 4, scale = TRUE)
PLSDA_3xtg_male_strain_scores <- data.frame(PLSDA_3xtg_male_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_male_strain_plot <- PLSDA_3xtg_male_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 3xTG Male - Strain",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_male_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_male_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_male_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_male_strain <- plotLoadings(PLSDA_3xtg_male_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_male_strain <- perf(PLSDA_3xtg_male_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_male_strain, legend = FALSE)

VIPs_3xtg_male_strain <- as.data.frame(mixOmics::vip(PLSDA_3xtg_male_strain))
VIPs_3xtg_male_strain_filter <- dplyr::filter(VIPs_3xtg_male_strain, VIPs_3xtg_male_strain$comp1 > 1)
VIPs_3xtg_male_strain_filter$ID <- rownames(VIPs_3xtg_male_strain_filter)
VIPs_3xtg_male_strain_select <- VIPs_3xtg_male_strain_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_male_strain_Load <- VIPs_3xtg_male_strain_select %>% 
  left_join(Loadings_3xtg_male_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 3xTG Male - Colonization
PLSDA_3xtg_male_colo <- mixOmics::plsda(data_3xtg_male_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                        PCA_3xtg_male_scores$Colonized, ncomp = 3, scale = TRUE)
PLSDA_3xtg_male_colo_scores <- data.frame(PLSDA_3xtg_male_colo$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_male_colo_plot <- PLSDA_3xtg_male_colo_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Colonized", alpha = 0.6, title = "PLSDA - Liver 3xTG Male - Colonized",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_male_colo$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_male_colo$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_male_colo_scores %>% group_by(Colonized) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Colonized), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_male_colo <- plotLoadings(PLSDA_3xtg_male_colo, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_male_colo <- perf(PLSDA_3xtg_male_colo, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_male_colo, legend = FALSE)

VIPs_3xtg_male_colo <- as.data.frame(mixOmics::vip(PLSDA_3xtg_male_colo))
VIPs_3xtg_male_colo_filter <- dplyr::filter(VIPs_3xtg_male_colo, VIPs_3xtg_male_colo$comp1 > 1)
VIPs_3xtg_male_colo_filter$ID <- rownames(VIPs_3xtg_male_colo_filter)
VIPs_3xtg_male_colo_select <- VIPs_3xtg_male_colo_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_male_colo_Load <- VIPs_3xtg_male_colo_select %>% 
  left_join(Loadings_3xtg_male_colo, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 3xTG Female - Strain
PLSDA_3xtg_female_strain <- mixOmics::plsda(data_3xtg_female_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                            PCA_3xtg_female_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_3xtg_female_strain_scores <- data.frame(PLSDA_3xtg_female_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_female_strain_plot <- PLSDA_3xtg_female_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 3xTG Female - Strain",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_female_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_female_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_female_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_female_strain <- plotLoadings(PLSDA_3xtg_female_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_female_strain <- perf(PLSDA_3xtg_female_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_female_strain, legend = FALSE)

VIPs_3xtg_female_strain <- as.data.frame(mixOmics::vip(PLSDA_3xtg_female_strain))
VIPs_3xtg_female_strain_filter <- dplyr::filter(VIPs_3xtg_female_strain, VIPs_3xtg_female_strain$comp1 > 1)
VIPs_3xtg_female_strain_filter$ID <- rownames(VIPs_3xtg_female_strain_filter)
VIPs_3xtg_female_strain_select <- VIPs_3xtg_female_strain_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_female_strain_Load <- VIPs_3xtg_female_strain_select %>% 
  left_join(Loadings_3xtg_female_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 3xTG Female - Colonization
PLSDA_3xtg_female_colo <- mixOmics::plsda(data_3xtg_female_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                          PCA_3xtg_female_scores$Colonized, ncomp = 3, scale = TRUE)
PLSDA_3xtg_female_colo_scores <- data.frame(PLSDA_3xtg_female_colo$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_female_colo_plot <- PLSDA_3xtg_female_colo_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Colonized", alpha = 0.6, title = "PLSDA - Liver 3xTG Female - Colonized",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_female_colo$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_female_colo$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_female_colo_scores %>% group_by(Colonized) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Colonized), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_female_colo <- plotLoadings(PLSDA_3xtg_female_colo, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_female_colo <- perf(PLSDA_3xtg_female_colo, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_female_colo, legend = FALSE)

VIPs_3xtg_female_colo <- as.data.frame(mixOmics::vip(PLSDA_3xtg_female_colo))
VIPs_3xtg_female_colo_filter <- dplyr::filter(VIPs_3xtg_female_colo, VIPs_3xtg_female_colo$comp1 > 1)
VIPs_3xtg_female_colo_filter$ID <- rownames(VIPs_3xtg_female_colo_filter)
VIPs_3xtg_female_colo_select <- VIPs_3xtg_female_colo_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_female_colo_Load <- VIPs_3xtg_female_colo_select %>% 
  left_join(Loadings_3xtg_female_colo, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# Sample matching 
list_3xtg_sex <- list(
  `Male 3xTG` = (VIPs_3xtg_male_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `Male WT` =  (VIPs_3xtg_male_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID,
  `Female 3xTG` = (VIPs_3xtg_female_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `Female WT` = (VIPs_3xtg_female_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID,
  `Male GF` = (VIPs_3xtg_male_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `Male SPF` =  (VIPs_3xtg_male_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID,
  `Female GF` = (VIPs_3xtg_female_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `Female SPF` = (VIPs_3xtg_female_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID)

sex_3xtg <- UpSetR::upset(fromList(list_3xtg_sex), nsets = 8, nintersects = NA, 
                          point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE, 
                          sets = c("Male 3xTG", "Male WT", "Female 3xTG", "Female WT",
                                   "Male GF", "Male SPF", "Female GF", "Female SPF"),
                          queries = list(list(query = intersects, params = list("Female SPF", 
                                                                                "Male SPF"), color = "steelblue1", active = T), 
                                         list(query = intersects, params = list("Female GF", 
                                                                                "Male GF"),color = "steelblue1", active = T), 
                                         list(query = intersects, params = list("Female 3xTG", 
                                                                                "Male 3xTG"),color = "orange", active = T), 
                                         list(query = intersects, params = list("Female WT", 
                                                                                "Male WT"),color = "orange", active = T)))


###################
# Combine microbiome and strain effect models to see intersection
###################
list_3xtg_intersection <- list(
  `3xTG - GF` = (VIPs_3xtg_mut_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `3xTG - SPF` = (VIPs_3xtg_mut_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID,
  `WT - GF` =  (VIPs_3xtg_wt_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `WT - SPF` = (VIPs_3xtg_wt_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID,
  `GF - 3xTG` = (VIPs_3xtg_gf_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `GF - WT` =  (VIPs_3xtg_gf_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID,
  `SPF - 3xTG` = (VIPs_3xtg_spf_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `SPF - WT` = (VIPs_3xtg_spf_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID)

intersection_3xtg <- UpSetR::upset(fromList(list_3xtg_intersection), nsets = 8, nintersects = NA, 
                                   point.size = 1.5, line.size = 1, text.scale = 1,
                                   keep.order = TRUE, 
                                   sets = c("3xTG - GF", "3xTG - SPF", "WT - GF", "WT - SPF",
                                            "GF - 3xTG", "GF - WT", "SPF - 3xTG", "SPF - WT"),
                                   queries = list(list(query = intersects, params = list("SPF - WT", 
                                                                                         "WT - SPF"), color = "steelblue1", active = T),
                                                  list(query = intersects, params = list("SPF - 3xTG", 
                                                                                         "3xTG - SPF"),color = "steelblue1", active = T),
                                                  list(query = intersects, params = list("SPF - WT", 
                                                                                         "GF - WT"),color = "green3", active = T),
                                                  list(query = intersects, params = list("SPF - 3xTG", 
                                                                                         "GF - 3xTG"),color = "green3", active = T),
                                                  list(query = intersects, params = list("WT - SPF", 
                                                                                         "3xTG - SPF"),color = "green4", active = T),
                                                  list(query = intersects, params = list("WT - GF", 
                                                                                         "3xTG - GF"),color = "green4", active = T),
                                                  list(query = intersects, params = list("GF - 3xTG", 
                                                                                         "3xTG - GF"),color = "orange", active = T),
                                                  list(query = intersects, params = list("GF - WT", 
                                                                                         "WT - GF"),color = "orange", active = T)))

#pdf(file="3xtg_liver_upset_genotype_colonization.pdf", width = 6.5, height = 4.5)
#intersection_3xtg
#dev.off()


# Extract IDs only present in `3xTG - SPF`, `SPF - 3xTG`, and `3xTG - SPF` U `SPF - 3xTG`
unique_3xtg_spf1 <- data.frame(ID  = setdiff(
  list_3xtg_intersection$`3xTG - SPF`,
  c(list_3xtg_intersection$`3xTG - GF`, list_3xtg_intersection$`WT - GF`, list_3xtg_intersection$`WT - SPF`,
    list_3xtg_intersection$`GF - 3xTG`, list_3xtg_intersection$`GF - WT`, list_3xtg_intersection$`SPF - 3xTG`, 
    list_3xtg_intersection$`SPF - WT`)), Group = "3xtg_spf")

unique_spf_3xtg1 <- data.frame(ID  = setdiff(
  list_3xtg_intersection$`SPF - 3xTG`,
  c(list_3xtg_intersection$`3xTG - GF`, list_3xtg_intersection$`WT - GF`, list_3xtg_intersection$`WT - SPF`,
    list_3xtg_intersection$`GF - 3xTG`, list_3xtg_intersection$`GF - WT`, list_3xtg_intersection$`3xTG - SPF`, 
    list_3xtg_intersection$`SPF - WT`)), Group = "spf_3xtg")

intersect_3xtg_spf_only <- intersect(
  list_3xtg_intersection$`3xTG - SPF`,
  list_3xtg_intersection$`SPF - 3xTG`)

unique_3xtg_spf_comb <- data.frame(ID  = setdiff(
  intersect_3xtg_spf_only,
  c(list_3xtg_intersection$`3xTG - GF`, list_3xtg_intersection$`WT - GF`, list_3xtg_intersection$`GF - 3xTG`, 
    list_3xtg_intersection$`GF - WT`, list_3xtg_intersection$`WT - SPF`, list_3xtg_intersection$`SPF - WT`)), Group = "3xtg_both")

combined_unique_3xtg_3xtg_spf <- rbind(unique_3xtg_spf1, unique_spf_3xtg1, unique_3xtg_spf_comb) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))

# Extract IDs only present in `WT - SPF`, `SPF - WT`, and `WT - SPF` U `SPF - WT`
unique_wt_spf1 <- data.frame(ID  = setdiff(
  list_3xtg_intersection$`WT - SPF`,
  c(list_3xtg_intersection$`3xTG - GF`, list_3xtg_intersection$`WT - GF`, list_3xtg_intersection$`3xTG - SPF`,
    list_3xtg_intersection$`GF - 3xTG`, list_3xtg_intersection$`GF - WT`, list_3xtg_intersection$`SPF - 3xTG`, 
    list_3xtg_intersection$`SPF - WT`)), Group = "wt_spf")

unique_spf_wt1 <- data.frame(ID = setdiff(
  list_3xtg_intersection$`SPF - WT`,
  c(list_3xtg_intersection$`3xTG - GF`, list_3xtg_intersection$`WT - GF`, list_3xtg_intersection$`WT - SPF`,
    list_3xtg_intersection$`GF - 3xTG`, list_3xtg_intersection$`GF - WT`, list_3xtg_intersection$`3xTG - SPF`, 
    list_3xtg_intersection$`SPF - 3xTG`)), Group = "spf_wt")

intersect_wt_spf_only <- intersect(
  list_3xtg_intersection$`WT - SPF`,
  list_3xtg_intersection$`SPF - WT`)

unique_wt_spf_comb <- data.frame(ID = setdiff(
  intersect_wt_spf_only,
  c(list_3xtg_intersection$`3xTG - GF`, list_3xtg_intersection$`WT - GF`, list_3xtg_intersection$`GF - 3xTG`, 
    list_3xtg_intersection$`GF - WT`, list_3xtg_intersection$`3xTG - SPF`, list_3xtg_intersection$`SPF - 3xTG`)), Group = "wt_both")

combined_unique_3xtg_wt_spf <- rbind(unique_wt_spf1, unique_spf_wt1, unique_wt_spf_comb) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))


# Extract feature of interest and check difference between strains colonized
exclusive_3xtg_liver <- rbind(combined_unique_3xtg_3xtg_spf, combined_unique_3xtg_wt_spf)

data_3xtg_feat <- data_sample %>% dplyr::filter(SampleID %in% sample_3xtg$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L3xtg_43", "L3xtg_55", "L3xtg_57", "L3xtg_59", 
                                  "L3xtg_49", "L3xtg_144", "L3xtg_117"))) %>%
  dplyr::select("SampleID", combined_unique_3xtg_3xtg_spf$ID, combined_unique_3xtg_wt_spf$ID) %>%
  dplyr::mutate(WT = rowSums(select(., combined_unique_3xtg_wt_spf$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., combined_unique_3xtg_3xtg_spf$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_metabolomics)

vip_3xtg_spf_wt <- VIPs_3xtg_spf_strain_Load %>% dplyr::filter(GroupContrib == "WT")
vip_3xtg_spf_mut <- VIPs_3xtg_spf_strain_Load %>% dplyr::filter(GroupContrib == "Mut")

data_3xtg_vip_spf <- data_sample %>% dplyr::filter(SampleID %in% sample_3xtg$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L3xtg_43", "L3xtg_55", "L3xtg_57", "L3xtg_59", 
                                  "L3xtg_49", "L3xtg_144", "L3xtg_117"))) %>%
  dplyr::select("SampleID", VIPs_3xtg_spf_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_3xtg_spf_wt$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_3xtg_spf_mut$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_metabolomics)

# Plot ratio
data_3xtg_feat %>% 
  dplyr::filter(Colonized == "SPF") %>%
  ggboxplot(y = "Ratio", x = "Strain", add = "jitter") + 
  stat_compare_means()

plot_ratio_liver_3xtg <- data_3xtg_vip_spf %>% 
  dplyr::filter(Colonized == "SPF") %>% 
  dplyr::mutate(Strain = factor(Strain, levels = c("WT", "Mut"))) %>%
  ggboxplot(y = "Ratio", x = "Strain", add = "jitter", ylab = "Ln(WT/Mut)",
            add.params = list(color = "Strain", alpha = 0.6), legend = "none",
            palette = c("#E69A8DFF", "#5F4B8BFF"), xlab = "Liver",
            title = "Differential features from PLS-DA models") + 
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = plot_ratio_liver_3xtg, filename = "3xtg_liver_ratio.svg", device = "svg", dpi = "retina", width = 1.5, height = 2.5)
#write_csv(x = data_3xtg_vip_spf, file = "liver_ratio_3xtg.csv")

# Cliff's Delta
cliff.delta((data_3xtg_vip_spf %>%
               dplyr::filter(Colonized == "SPF" & Strain == "WT"))$Ratio, 
            (data_3xtg_vip_spf %>%
               dplyr::filter(Colonized == "SPF" & Strain == "Mut"))$Ratio)


# Extract features of interest for MN
liver_3xtg_interest <- VIPs_3xtg_spf_strain_Load %>% 
  dplyr::select(ID, comp1, GroupContrib) %>% 
  left_join(VIPs_3xtg_gf_strain_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>% 
  left_join(VIPs_3xtg_mut_colo_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>% 
  left_join(VIPs_3xtg_wt_colo_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>% 
  dplyr::mutate(VIP_mean = rowMeans(dplyr::select(., comp1.x, comp1.x.x, comp1.y, comp1.y.y), na.rm = TRUE)) %>%
  left_join(exclusive_3xtg_liver %>% dplyr::select(ID, Group)) %>% arrange(desc(`comp1.x`))

colnames(liver_3xtg_interest) <- c("ID", "VIP1", "SPF_genotype", "VIP2", "GF_genotype",
                                   "VIP3", "3xTG_colonization", "VIP4", "WT_colonization", 
                                   "VIP_mean", "Exclusivity")

liver_3xtg_interest <- liver_3xtg_interest %>%   
  dplyr::mutate(Genotype = case_when(SPF_genotype == GF_genotype ~ "Yes",
                                     TRUE ~ "No")) %>%
  dplyr::mutate(Microbiome = case_when(`3xTG_colonization` == WT_colonization ~ "Yes",
                                       TRUE ~ "No")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))

#write_csv(x = liver_3xtg_interest, file = "3xtg_liver_features_info_mn.csv")


#########
# 5xFAD #
#########
data_5xfad <- data_sample %>% dplyr::filter(SampleID %in% sample_5xfad$SampleID)

# RCLR transformation
data_5xfad_clr <- decostand(data_5xfad %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_5xfad_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

i <- "box"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Liver 5xFAD", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Some of the samples from Box 9 seems to be compromised. Remove them
samples_to_remove <- PCA_whole_scores %>% dplyr::filter(PC1 < 0)

data_5xfad_filter <- data_5xfad %>% dplyr::filter(!(SampleID %in% samples_to_remove$SampleID)) %>%
  dplyr::filter(!(SampleID %in% c("L5xfad_134", "L5xfad_98")))

# RCLR transformation
data_5xfad_clr <- decostand(data_5xfad_filter %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_5xfad_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PCA_5xfad_plots <- list()

for (i in c("box", "host_age", "sex", "Plate",  "Colonized", "Strain")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 5xFAD", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_5xfad_plots[[i]] <- PCA_plot
  
}

PCA_5xfad_plots_final <- wrap_plots(PCA_5xfad_plots, nrow = 3)

# PERMANOVA
dist_metabolites <- vegdist(data_5xfad_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$mouse_strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ Strain * Colonized * sex + 
                       box + Plate, PCA_whole_scores, na.action = na.omit)


# PLSDA 5xFAD - Strain
PLSDA_5xfad_strain <- mixOmics::plsda(data_5xfad_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                      PCA_whole_scores$Strain, ncomp = 2, scale = TRUE)
PLSDA_5xfad_strain_scores <- data.frame(PLSDA_5xfad_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_strain_plot <- PLSDA_5xfad_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 5xFAD - Strain",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_strain <- plotLoadings(PLSDA_5xfad_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_strain <- perf(PLSDA_5xfad_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_strain, legend = FALSE)

VIPs_5xfad_strain <- as.data.frame(mixOmics::vip(PLSDA_5xfad_strain))
VIPs_5xfad_strain_filter <- dplyr::filter(VIPs_5xfad_strain, VIPs_5xfad_strain$comp1 > 1)
VIPs_5xfad_strain_filter$ID <- rownames(VIPs_5xfad_strain_filter)
VIPs_5xfad_strain_select <- VIPs_5xfad_strain_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_strain_Load <- VIPs_5xfad_strain_select %>% 
  left_join(Loadings_5xfad_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 5xFAD - Colonized
PLSDA_5xfad_colonized <- mixOmics::plsda(data_5xfad_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                         PCA_whole_scores$Colonized, ncomp = 2, scale = TRUE)
PLSDA_5xfad_colonized_scores <- data.frame(PLSDA_5xfad_colonized$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_colonized_plot <- PLSDA_5xfad_colonized_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Colonized", alpha = 0.6, title = "PLSDA - Liver 5xFAD - Colonized",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_colonized$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_colonized$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_colonized_scores %>% group_by(Colonized) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Colonized), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_colonized <- plotLoadings(PLSDA_5xfad_colonized, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_colonized <- perf(PLSDA_5xfad_colonized, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_colonized, legend = FALSE)

VIPs_5xfad_colo <- as.data.frame(mixOmics::vip(PLSDA_5xfad_colonized))
VIPs_5xfad_colo_filter <- dplyr::filter(VIPs_5xfad_colo, VIPs_5xfad_colo$comp1 > 1)
VIPs_5xfad_colo_filter$ID <- rownames(VIPs_5xfad_colo_filter)
VIPs_5xfad_colo_select <- VIPs_5xfad_colo_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_colo_Load <- VIPs_5xfad_colo_select %>% 
  left_join(Loadings_5xfad_colonized, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 5xFAD - Sex
PLSDA_5xfad_sex <- mixOmics::plsda(data_5xfad_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                   PCA_whole_scores$sex, ncomp = 2, scale = TRUE)
PLSDA_5xfad_sex_scores <- data.frame(PLSDA_5xfad_sex$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_sex_plot <- PLSDA_5xfad_sex_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "sex", alpha = 0.6, title = "PLSDA - Liver 5xFAD - Sex",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_sex$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_sex$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_sex_scores %>% group_by(sex) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = sex), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_sex <- plotLoadings(PLSDA_5xfad_sex, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_sex <- perf(PLSDA_5xfad_sex, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_sex, legend = FALSE)

VIPs_5xfad_sex <- as.data.frame(mixOmics::vip(PLSDA_5xfad_sex))
VIPs_5xfad_sex_filter <- dplyr::filter(VIPs_5xfad_sex, VIPs_5xfad_sex$comp1 > 1)
VIPs_5xfad_sex_filter$ID <- rownames(VIPs_5xfad_sex_filter)
VIPs_5xfad_sex_select <- VIPs_5xfad_sex_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_sex_Load <- VIPs_5xfad_sex_select %>% 
  left_join(Loadings_5xfad_sex, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))


# Remove feature influenced by sex
VIPs_5xfad_strain_no_sex <- VIPs_5xfad_strain_Load %>% 
  anti_join(VIPs_5xfad_sex_Load %>% dplyr::select(ID))
VIPs_5xfad_colo_no_sex <- VIPs_5xfad_colo_Load %>% 
  anti_join(VIPs_5xfad_sex_Load %>% dplyr::select(ID))

VIPs_5xfad_colo_strain_clean <- VIPs_5xfad_strain_no_sex %>% 
  inner_join(VIPs_5xfad_colo_no_sex %>% dplyr::select(ID))

#write_csv(x = VIPs_5xfad_strain_no_sex, file = "5xfad_liver_plsda_vip_all_samples_genotype.csv")
#write_csv(x = VIPs_5xfad_colo_no_sex, file = "5xfad_liver_plsda_vip_all_samples_colonization.csv")
#write_csv(x = VIPs_5xfad_colo_strain_clean, file = "5xfad_liver_plsda_vip_all_samples_genotype_colonization_intersection.csv")


############
# Compare 5xfad gf v wt gf to 5xfad spf v wt spf to see genotype effect
############
data_5xfad_gf <- data_5xfad %>% dplyr::filter(!(SampleID %in% samples_to_remove$SampleID)) %>%
  dplyr::filter(SampleID %in% sample_gf$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L5xfad_134", "L5xfad_98")))

data_5xfad_spf <- data_5xfad %>% dplyr::filter(!(SampleID %in% samples_to_remove$SampleID)) %>%
  dplyr::filter(SampleID %in% sample_spf$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L5xfad_134", "L5xfad_98")))

# RCLR transformation
data_5xfad_gf_clr <- decostand(data_5xfad_gf %>% column_to_rownames("SampleID"), method = "rclr")
data_5xfad_spf_clr <- decostand(data_5xfad_spf %>% column_to_rownames("SampleID"), method = "rclr")

# PCA GF
PCA_5xfad_gf <- mixOmics::pca(data_5xfad_gf_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                              ncomp = 2, center = TRUE, scale = TRUE)
PCA_5xfad_gf_scores <- data.frame(PCA_5xfad_gf$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_5xfad_gf_plots <- list()

for (i in c("box", "host_age", "sex", "Plate", "Strain")) {
  
  PCA_plot <- PCA_5xfad_gf_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 5xFAD - GF", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_5xfad_gf$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_5xfad_gf$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_5xfad_gf_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_5xfad_gf_plots[[i]] <- PCA_plot
  
}

PCA_5xfad_gf_plots_final <- wrap_plots(PCA_5xfad_gf_plots, nrow = 3)

# PCA SPF
PCA_5xfad_spf <- mixOmics::pca(data_5xfad_spf_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                               ncomp = 2, center = TRUE, scale = TRUE)
PCA_5xfad_spf_scores <- data.frame(PCA_5xfad_spf$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_5xfad_spf_plots <- list()

for (i in c("box", "host_age", "sex", "Plate", "Strain")) {
  
  PCA_plot <- PCA_5xfad_spf_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 5xFAD SPF", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_5xfad_spf$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_5xfad_spf$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_5xfad_spf_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_5xfad_spf_plots[[i]] <- PCA_plot
  
}

PCA_5xfad_spf_plots_final <- wrap_plots(PCA_5xfad_spf_plots, nrow = 3)


# PLSDA 5xFAD GF - Strain
PLSDA_5xfad_gf_strain <- mixOmics::plsda(data_5xfad_gf_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                         PCA_5xfad_gf_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_5xfad_gf_strain_scores <- data.frame(PLSDA_5xfad_gf_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_gf_strain_plot <- PLSDA_5xfad_gf_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 5xFAD GF - Strain",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_gf_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_gf_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_gf_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_gf_strain <- plotLoadings(PLSDA_5xfad_gf_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_gf_strain <- perf(PLSDA_5xfad_gf_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_gf_strain, legend = FALSE)

VIPs_5xfad_gf_strain <- as.data.frame(mixOmics::vip(PLSDA_5xfad_gf_strain))
VIPs_5xfad_gf_strain_filter <- dplyr::filter(VIPs_5xfad_gf_strain, VIPs_5xfad_gf_strain$comp1 > 1)
VIPs_5xfad_gf_strain_filter$ID <- rownames(VIPs_5xfad_gf_strain_filter)
VIPs_5xfad_gf_strain_select <- VIPs_5xfad_gf_strain_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_gf_strain_Load <- VIPs_5xfad_gf_strain_select %>% 
  left_join(Loadings_5xfad_gf_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 5xFAD SPF - Strain
PLSDA_5xfad_spf_strain <- mixOmics::plsda(data_5xfad_spf_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                          PCA_5xfad_spf_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_5xfad_spf_strain_scores <- data.frame(PLSDA_5xfad_spf_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_spf_strain_plot <- PLSDA_5xfad_spf_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 5xFAD SPF - Strain",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_spf_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_spf_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_spf_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_spf_strain <- plotLoadings(PLSDA_5xfad_spf_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_spf_strain <- perf(PLSDA_5xfad_spf_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_spf_strain, legend = FALSE)

#write_csv(x = perf_plsda_5xfad_spf_strain$error.rate$overall %>% 
#            as.data.frame() %>% dplyr::mutate(Model = "5xfad_spf_strain"), file = "liver_plsda_performance_5xfad.csv")

VIPs_5xfad_spf_strain <- as.data.frame(mixOmics::vip(PLSDA_5xfad_spf_strain))
VIPs_5xfad_spf_strain_filter <- dplyr::filter(VIPs_5xfad_spf_strain, VIPs_5xfad_spf_strain$comp1 > 1)
VIPs_5xfad_spf_strain_filter$ID <- rownames(VIPs_5xfad_spf_strain_filter)
VIPs_5xfad_spf_strain_select <- VIPs_5xfad_spf_strain_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_spf_strain_Load <- VIPs_5xfad_spf_strain_select %>% 
  left_join(Loadings_5xfad_spf_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# Sample matching 
list_5xfad_strain <- list(
  `GF 5xFAD` = (VIPs_5xfad_gf_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `GF WT` =  (VIPs_5xfad_gf_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID,
  `SPF 5xFAD` = (VIPs_5xfad_spf_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `SPF WT` = (VIPs_5xfad_spf_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID)

strain_5xfad <- UpSetR::upset(fromList(list_5xfad_strain), nsets = 4, nintersects = 6, 
                              point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                              sets = c("GF 5xFAD", "GF WT", "SPF 5xFAD", "SPF WT"),
                              queries = list(list(query = intersects, params = list("SPF 5xFAD", 
                                                                                    "GF 5xFAD"), color = "steelblue1", active = T), 
                                             list(query = intersects, params = list("SPF WT", 
                                                                                    "GF WT"),color = "orange", active = T)))

#pdf(file="5xfad_liver_upset_genotype.pdf", width = 3, height = 2.5)
#strain_5xfad
#dev.off()


# Check common features 5xfad mutation or wt
features_5xfad_mut <- VIPs_5xfad_spf_strain_Load %>% 
  dplyr::filter(GroupContrib == "Mut") %>%
  dplyr::select(ID, comp1) %>% 
  inner_join(VIPs_5xfad_gf_strain_Load %>% 
               dplyr::filter(GroupContrib == "Mut"), by = c("ID" = "ID"))

features_5xfad_wt <- VIPs_5xfad_spf_strain_Load %>% 
  dplyr::filter(GroupContrib == "WT") %>%
  dplyr::select(ID, comp1) %>% 
  inner_join(VIPs_5xfad_gf_strain_Load %>% 
               dplyr::filter(GroupContrib == "WT"), by = c("ID" = "ID"))

#write_csv(x = features_5xfad_mut, file = "5xfad_liver_genotype_mut.csv")
#write_csv(x = features_5xfad_wt, file = "5xfad_liver_genotype_wt.csv")


############
# Compare 5xfad gf v 5xfad spf to wt gf v wt spf to see microbiome effect on strain
############
data_5xfad_mut <- data_5xfad %>% dplyr::filter(!(SampleID %in% samples_to_remove$SampleID)) %>%
  dplyr::filter(SampleID %in% sample_mut$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L5xfad_134", "L5xfad_98")))

data_5xfad_wt <- data_5xfad %>% dplyr::filter(!(SampleID %in% samples_to_remove$SampleID)) %>%
  dplyr::filter(SampleID %in% sample_wt$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L5xfad_134", "L5xfad_98")))

# RCLR transformation
data_5xfad_mut_clr <- decostand(data_5xfad_mut %>% column_to_rownames("SampleID"), method = "rclr")
data_5xfad_wt_clr <- decostand(data_5xfad_wt %>% column_to_rownames("SampleID"), method = "rclr")

# PCA Mut
PCA_5xfad_mut <- mixOmics::pca(data_5xfad_mut_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                               ncomp = 2, center = TRUE, scale = TRUE)
PCA_5xfad_mut_scores <- data.frame(PCA_5xfad_mut$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_5xfad_mut_plots <- list()

for (i in c("box", "host_age", "Plate", "Colonized", "sex")) {
  
  PCA_plot <- PCA_5xfad_mut_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 5xFAD - Mut", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_5xfad_mut$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_5xfad_mut$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_5xfad_mut_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_5xfad_mut_plots[[i]] <- PCA_plot
  
}

PCA_5xfad_mut_plots_final <- wrap_plots(PCA_5xfad_mut_plots, nrow = 3)

# PCA WT
PCA_5xfad_wt <- mixOmics::pca(data_5xfad_wt_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                              ncomp = 2, center = TRUE, scale = TRUE)
PCA_5xfad_wt_scores <- data.frame(PCA_5xfad_wt$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_5xfad_wt_plots <- list()

for (i in c("box", "host_age", "Plate", "Colonized", "sex")) {
  
  PCA_plot <- PCA_5xfad_wt_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 5xFAD - WT", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_5xfad_wt$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_5xfad_wt$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_5xfad_wt_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_5xfad_wt_plots[[i]] <- PCA_plot
  
}

PCA_5xfad_wt_plots_final <- wrap_plots(PCA_5xfad_wt_plots, nrow = 3)


# PLSDA 5xFAD Mut - Colonization
PLSDA_5xfad_mut_colo <- mixOmics::plsda(data_5xfad_mut_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                          PCA_5xfad_mut_scores$Colonized, ncomp = 2, scale = TRUE)
PLSDA_5xfad_mut_colo_scores <- data.frame(PLSDA_5xfad_mut_colo$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_mut_colo_plot <- PLSDA_5xfad_mut_colo_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Colonized", alpha = 0.6, title = "PLSDA - Liver 5xFAD Mut - Colonization",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_mut_colo$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_mut_colo$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_mut_colo_scores %>% group_by(Colonized) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Colonized), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_mut_strain <- plotLoadings(PLSDA_5xfad_mut_colo, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_mut_colonization <- perf(PLSDA_5xfad_mut_colo, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_mut_colonization, legend = FALSE)

VIPs_5xfad_mut_colo <- as.data.frame(mixOmics::vip(PLSDA_5xfad_mut_colo))
VIPs_5xfad_mut_colo_filter <- dplyr::filter(VIPs_5xfad_mut_colo, VIPs_5xfad_mut_colo$comp1 > 1)
VIPs_5xfad_mut_colo_filter$ID <- rownames(VIPs_5xfad_mut_colo_filter)
VIPs_5xfad_mut_colo_select <- VIPs_5xfad_mut_colo_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_mut_colo_Load <- VIPs_5xfad_mut_colo_select %>% 
  left_join(Loadings_5xfad_mut_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 5xFAD WT - Colonization
PLSDA_5xfad_wt_colo <- mixOmics::plsda(data_5xfad_wt_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                         PCA_5xfad_wt_scores$Colonized, ncomp = 2, scale = TRUE)
PLSDA_5xfad_wt_colo_scores <- data.frame(PLSDA_5xfad_wt_colo$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_wt_colo_plot <- PLSDA_5xfad_wt_colo_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Colonized", alpha = 0.6, title = "PLSDA - Liver 5xFAD WT - Colonization",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_wt_colo$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_wt_colo$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_wt_colo_scores %>% group_by(Colonized) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Colonized), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_wt_strain <- plotLoadings(PLSDA_5xfad_wt_colo, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_wt_colonization <- perf(PLSDA_5xfad_wt_colo, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_wt_colonization, legend = FALSE)

VIPs_5xfad_wt_colo <- as.data.frame(mixOmics::vip(PLSDA_5xfad_wt_colo))
VIPs_5xfad_wt_colo_filter <- dplyr::filter(VIPs_5xfad_wt_colo, VIPs_5xfad_wt_colo$comp1 > 1)
VIPs_5xfad_wt_colo_filter$ID <- rownames(VIPs_5xfad_wt_colo_filter)
VIPs_5xfad_wt_colo_select <- VIPs_5xfad_wt_colo_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_wt_colo_Load <- VIPs_5xfad_wt_colo_select %>% 
  left_join(Loadings_5xfad_wt_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# Sample matching 
list_5xfad_microbiome <- list(
  `5xFAD GF` = (VIPs_5xfad_mut_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `5xFAD SPF` = (VIPs_5xfad_mut_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID,
  `WT GF` =  (VIPs_5xfad_wt_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `WT SPF` = (VIPs_5xfad_wt_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID)

colonization_5xfad <- UpSetR::upset(fromList(list_5xfad_microbiome), nsets = 4, nintersects = 6, 
                                    point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                                    sets = c("5xFAD GF", "5xFAD SPF", "WT GF", "WT SPF"),
                                    queries = list(list(query = intersects, params = list("WT SPF", 
                                                                                          "5xFAD SPF"), color = "steelblue1", active = T), 
                                                   list(query = intersects, params = list("WT GF", 
                                                                                          "5xFAD GF"),color = "orange", active = T)))

#pdf(file="5xfad_liver_upset_colonization.pdf", width = 3, height = 2.5)
#colonization_5xfad
#dev.off()

# Check common features spf or gf
features_5xfad_spf <- VIPs_5xfad_mut_colo_Load %>% 
  dplyr::filter(GroupContrib == "SPF") %>%
  dplyr::select(ID, comp1) %>% 
  inner_join(VIPs_5xfad_wt_colo_Load %>% 
               dplyr::filter(GroupContrib == "SPF"), by = c("ID" = "ID"))

features_5xfad_gf <- VIPs_5xfad_mut_colo_Load %>% 
  dplyr::filter(GroupContrib == "GF") %>%
  dplyr::select(ID, comp1) %>% 
  inner_join(VIPs_5xfad_wt_colo_Load %>% 
               dplyr::filter(GroupContrib == "GF"), by = c("ID" = "ID"))

#write_csv(x = features_5xfad_spf, file = "5xfad_liver_colonization_spf.csv")
#write_csv(x = features_5xfad_gf, file = "5xfad_liver_colonization_gf.csv")


############
# Compare 5xfad male v wt male to 5xfad female v wt female to see strain effect on sex
############
data_5xfad_male <- data_5xfad %>% dplyr::filter(!(SampleID %in% samples_to_remove$SampleID)) %>%
  dplyr::filter(SampleID %in% sample_male$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L5xfad_134", "L5xfad_98")))

data_5xfad_female <- data_5xfad %>% dplyr::filter(!(SampleID %in% samples_to_remove$SampleID)) %>%
  dplyr::filter(SampleID %in% sample_female$SampleID) %>%
  dplyr::filter(!(SampleID %in% c("L5xfad_134", "L5xfad_98")))

# RCLR transformation
data_5xfad_male_clr <- decostand(data_5xfad_male %>% column_to_rownames("SampleID"), method = "rclr")
data_5xfad_female_clr <- decostand(data_5xfad_female %>% column_to_rownames("SampleID"), method = "rclr")

# PCA Male
PCA_5xfad_male <- mixOmics::pca(data_5xfad_male_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                ncomp = 2, center = TRUE, scale = TRUE)
PCA_5xfad_male_scores <- data.frame(PCA_5xfad_male$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_5xfad_male_plots <- list()

for (i in c("box", "host_age", "Plate", "Strain", "Colonized")) {
  
  PCA_plot <- PCA_5xfad_male_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 5xFAD Male", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_5xfad_male$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_5xfad_male$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_5xfad_male_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_5xfad_male_plots[[i]] <- PCA_plot
  
}

PCA_5xfad_male_plots_final <- wrap_plots(PCA_5xfad_male_plots, nrow = 3)

# PCA Female
PCA_5xfad_female <- mixOmics::pca(data_5xfad_female_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                  ncomp = 2, center = TRUE, scale = TRUE)
PCA_5xfad_female_scores <- data.frame(PCA_5xfad_female$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_5xfad_female_plots <- list()

for (i in c("box", "host_age", "Plate", "Strain", "Colonized")) {
  
  PCA_plot <- PCA_5xfad_female_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Liver 5xFAD Female", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_5xfad_female$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_5xfad_female$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_5xfad_female_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_5xfad_female_plots[[i]] <- PCA_plot
  
}

PCA_5xfad_female_plots_final <- wrap_plots(PCA_5xfad_female_plots, nrow = 3)


# PLSDA 5xFAD Male - Strain
PLSDA_5xfad_male_strain <- mixOmics::plsda(data_5xfad_male_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                           PCA_5xfad_male_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_5xfad_male_strain_scores <- data.frame(PLSDA_5xfad_male_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_male_strain_plot <- PLSDA_5xfad_male_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 5xFAD Male - Strain",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_male_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_male_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_male_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_male_strain <- plotLoadings(PLSDA_5xfad_male_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_male_strain <- perf(PLSDA_5xfad_male_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_male_strain, legend = FALSE)

VIPs_5xfad_male_strain <- as.data.frame(mixOmics::vip(PLSDA_5xfad_male_strain))
VIPs_5xfad_male_strain_filter <- dplyr::filter(VIPs_5xfad_male_strain, VIPs_5xfad_male_strain$comp1 > 1)
VIPs_5xfad_male_strain_filter$ID <- rownames(VIPs_5xfad_male_strain_filter)
VIPs_5xfad_male_strain_select <- VIPs_5xfad_male_strain_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_male_strain_Load <- VIPs_5xfad_male_strain_select %>% 
  left_join(Loadings_5xfad_male_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 5xFAD Male - Colonization
PLSDA_5xfad_male_colo <- mixOmics::plsda(data_5xfad_male_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                         PCA_5xfad_male_scores$Colonized, ncomp = 3, scale = TRUE)
PLSDA_5xfad_male_colo_scores <- data.frame(PLSDA_5xfad_male_colo$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_male_colo_plot <- PLSDA_5xfad_male_colo_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Colonized", alpha = 0.6, title = "PLSDA - Liver 5xFAD Male - Colonized",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_male_colo$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_male_colo$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_male_colo_scores %>% group_by(Colonized) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Colonized), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_male_colo <- plotLoadings(PLSDA_5xfad_male_colo, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_male_colo <- perf(PLSDA_5xfad_male_colo, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_male_colo, legend = FALSE)

VIPs_5xfad_male_colo <- as.data.frame(mixOmics::vip(PLSDA_5xfad_male_colo))
VIPs_5xfad_male_colo_filter <- dplyr::filter(VIPs_5xfad_male_colo, VIPs_5xfad_male_colo$comp1 > 1)
VIPs_5xfad_male_colo_filter$ID <- rownames(VIPs_5xfad_male_colo_filter)
VIPs_5xfad_male_colo_select <- VIPs_5xfad_male_colo_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_male_colo_Load <- VIPs_5xfad_male_colo_select %>% 
  left_join(Loadings_5xfad_male_colo, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 5xFAD Female - Strain
PLSDA_5xfad_female_strain <- mixOmics::plsda(data_5xfad_female_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                             PCA_5xfad_female_scores$Strain, ncomp = 2, scale = TRUE)
PLSDA_5xfad_female_strain_scores <- data.frame(PLSDA_5xfad_female_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_female_strain_plot <- PLSDA_5xfad_female_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 5xFAD Female - Strain",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_female_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_female_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_female_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_female_strain <- plotLoadings(PLSDA_5xfad_female_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_female_strain <- perf(PLSDA_5xfad_female_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_female_strain, legend = FALSE)

VIPs_5xfad_female_strain <- as.data.frame(mixOmics::vip(PLSDA_5xfad_female_strain))
VIPs_5xfad_female_strain_filter <- dplyr::filter(VIPs_5xfad_female_strain, VIPs_5xfad_female_strain$comp1 > 1)
VIPs_5xfad_female_strain_filter$ID <- rownames(VIPs_5xfad_female_strain_filter)
VIPs_5xfad_female_strain_select <- VIPs_5xfad_female_strain_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_female_strain_Load <- VIPs_5xfad_female_strain_select %>% 
  left_join(Loadings_5xfad_female_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 5xFAD Female - Colonization
PLSDA_5xfad_female_colo <- mixOmics::plsda(data_5xfad_female_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                           PCA_5xfad_female_scores$Colonized, ncomp = 3, scale = TRUE)
PLSDA_5xfad_female_colo_scores <- data.frame(PLSDA_5xfad_female_colo$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_female_colo_plot <- PLSDA_5xfad_female_colo_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Liver 5xFAD Female - Colonized",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_female_colo$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_female_colo$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_female_colo_scores %>% group_by(Colonized) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Colonized), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_female_colo <- plotLoadings(PLSDA_5xfad_female_colo, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_female_colo <- perf(PLSDA_5xfad_female_colo, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_female_colo, legend = FALSE)

VIPs_5xfad_female_colo <- as.data.frame(mixOmics::vip(PLSDA_5xfad_female_colo))
VIPs_5xfad_female_colo_filter <- dplyr::filter(VIPs_5xfad_female_colo, VIPs_5xfad_female_colo$comp1 > 1)
VIPs_5xfad_female_colo_filter$ID <- rownames(VIPs_5xfad_female_colo_filter)
VIPs_5xfad_female_colo_select <- VIPs_5xfad_female_colo_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_female_colo_Load <- VIPs_5xfad_female_colo_select %>% 
  left_join(Loadings_5xfad_female_colo, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# Sample matching 
list_5xfad_sex <- list(
  `Male 5xFAD` = (VIPs_5xfad_male_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `Male WT` =  (VIPs_5xfad_male_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID,
  `Female 5xFAD` = (VIPs_5xfad_female_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `Female WT` = (VIPs_5xfad_female_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID,
  `Male GF` = (VIPs_5xfad_male_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `Male SPF` =  (VIPs_5xfad_male_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID,
  `Female GF` = (VIPs_5xfad_female_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `Female SPF` = (VIPs_5xfad_female_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID)

sex_5xfad <- UpSetR::upset(fromList(list_5xfad_sex), nsets = 8, nintersects = NA, 
                           point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE, 
                           sets = c("Male 5xFAD", "Male WT", "Female 5xFAD", "Female WT",
                                    "Male GF", "Male SPF", "Female GF", "Female SPF"),
                           queries = list(list(query = intersects, params = list("Female SPF", 
                                                                                 "Male SPF"), color = "steelblue1", active = T), 
                                          list(query = intersects, params = list("Female GF", 
                                                                                 "Male GF"),color = "steelblue1", active = T), 
                                          list(query = intersects, params = list("Female 5xFAD", 
                                                                                 "Male 5xFAD"),color = "orange", active = T), 
                                          list(query = intersects, params = list("Female WT", 
                                                                                 "Male WT"),color = "orange", active = T)))


###################
# Combine microbiome and strain effect models to see intersection
###################
list_5xfad_intersection <- list(
  `5xFAD - GF` = (VIPs_5xfad_mut_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `5xFAD - SPF` = (VIPs_5xfad_mut_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID,
  `WT - GF` =  (VIPs_5xfad_wt_colo_Load %>% dplyr::filter(GroupContrib == "GF"))$ID,
  `WT - SPF` = (VIPs_5xfad_wt_colo_Load %>% dplyr::filter(GroupContrib == "SPF"))$ID,
  `GF - 5xFAD` = (VIPs_5xfad_gf_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `GF - WT` =  (VIPs_5xfad_gf_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID,
  `SPF - 5xFAD` = (VIPs_5xfad_spf_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `SPF - WT` = (VIPs_5xfad_spf_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID)

intersection_5xfad <- UpSetR::upset(fromList(list_5xfad_intersection), nsets = 8, nintersects = NA, 
                                    point.size = 1.5, line.size = 1, text.scale = 1,
                                    keep.order = TRUE, 
                                    sets = c("5xFAD - GF", "5xFAD - SPF", "WT - GF", "WT - SPF",
                                             "GF - 5xFAD", "GF - WT", "SPF - 5xFAD", "SPF - WT"),
                                    queries = list(list(query = intersects, params = list("SPF - WT", 
                                                                                          "WT - SPF"), color = "steelblue1", active = T),
                                                   list(query = intersects, params = list("SPF - 5xFAD", 
                                                                                          "5xFAD - SPF"),color = "steelblue1", active = T),
                                                   list(query = intersects, params = list("SPF - WT", 
                                                                                          "GF - WT"),color = "green3", active = T),
                                                   list(query = intersects, params = list("SPF - 5xFAD", 
                                                                                          "GF - 5xFAD"),color = "green3", active = T),
                                                   list(query = intersects, params = list("WT - SPF", 
                                                                                          "5xFAD - SPF"),color = "green4", active = T),
                                                   list(query = intersects, params = list("WT - GF", 
                                                                                          "5xFAD - GF"),color = "green4", active = T),
                                                   list(query = intersects, params = list("GF - 5xFAD", 
                                                                                          "5xFAD - GF"),color = "orange", active = T),
                                                   list(query = intersects, params = list("GF - WT", 
                                                                                          "WT - GF"),color = "orange", active = T)))


#pdf(file="5xfad_liver_upset_genotype_colonization.pdf", width = 6.5, height = 4.5)
#intersection_5xfad
#dev.off()

# Extract IDs only present in `5xFAD - SPF`, `SPF - 5xFAD`, and 5xFAD - SPF` U `SPF - 5xFAD`
unique_5xfad_spf1 <- data.frame(ID  = setdiff(
  list_5xfad_intersection$`5xFAD - SPF`,
  c(list_5xfad_intersection$`5xFAD - GF`, list_5xfad_intersection$`WT - GF`, list_5xfad_intersection$`WT - SPF`,
    list_5xfad_intersection$`GF - 5xFAD`, list_5xfad_intersection$`GF - WT`, list_5xfad_intersection$`SPF - 5xFAD`, 
    list_5xfad_intersection$`SPF - WT`)), Group = "5xfad_spf")

unique_spf_5xfad1 <- data.frame(ID  = setdiff(
  list_5xfad_intersection$`SPF - 5xFAD`,
  c(list_5xfad_intersection$`5xFAD - GF`, list_5xfad_intersection$`WT - GF`, list_5xfad_intersection$`WT - SPF`,
    list_5xfad_intersection$`GF - 5xFAD`, list_5xfad_intersection$`GF - WT`, list_5xfad_intersection$`5xFAD - SPF`, 
    list_5xfad_intersection$`SPF - WT`)), Group = "spf_5xfad")

intersect_5xfad_spf_only <- intersect(
  list_5xfad_intersection$`5xFAD - SPF`,
  list_5xfad_intersection$`SPF - 5xFAD`)

unique_5xfad_spf_comb <- data.frame(ID  = setdiff(
  intersect_5xfad_spf_only,
  c(list_5xfad_intersection$`5xFAD - GF`, list_5xfad_intersection$`WT - GF`, list_5xfad_intersection$`GF - 5xFAD`, 
    list_5xfad_intersection$`GF - WT`, list_5xfad_intersection$`WT - SPF`, list_5xfad_intersection$`SPF - WT`)), Group = "5xfad_both")

combined_unique_5xfad_5xfad_spf <- rbind(unique_5xfad_spf1, unique_spf_5xfad1, unique_5xfad_spf_comb) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))

# Extract IDs only present in `WT - SPF`, `SPF - WT`, and `WT - SPF` U `SPF - WT`
unique_wt_5xfad_spf1 <- data.frame(ID  = setdiff(
  list_5xfad_intersection$`WT - SPF`,
  c(list_5xfad_intersection$`5xFAD - GF`, list_5xfad_intersection$`WT - GF`, list_5xfad_intersection$`5xFAD - SPF`,
    list_5xfad_intersection$`GF - 5xFAD`, list_5xfad_intersection$`GF - WT`, list_5xfad_intersection$`SPF - 5xFAD`, 
    list_5xfad_intersection$`SPF - WT`)), Group = "wt_spf")

unique_spf_5xfad_wt1 <- data.frame(ID = setdiff(
  list_5xfad_intersection$`SPF - WT`,
  c(list_5xfad_intersection$`5xFAD - GF`, list_5xfad_intersection$`WT - GF`, list_5xfad_intersection$`WT - SPF`,
    list_5xfad_intersection$`GF - 5xFAD`, list_5xfad_intersection$`GF - WT`, list_5xfad_intersection$`5xFAD - SPF`, 
    list_5xfad_intersection$`SPF - 5xFAD`)), Group = "spf_wt")

intersect_wt_spf_5xfad_only <- intersect(
  list_5xfad_intersection$`WT - SPF`,
  list_5xfad_intersection$`SPF - WT`)

unique_wt_spf_5xfad_comb <- data.frame(ID = setdiff(
  intersect_wt_spf_5xfad_only,
  c(list_5xfad_intersection$`5xFAD - GF`, list_5xfad_intersection$`WT - GF`, list_5xfad_intersection$`GF - 5xFAD`, 
    list_5xfad_intersection$`GF - WT`, list_5xfad_intersection$`5xFAD - SPF`, list_5xfad_intersection$`SPF - 5xFAD`)), Group = "wt_both")

combined_unique_5xfad_wt_spf <- rbind(unique_wt_5xfad_spf1, unique_spf_5xfad_wt1, unique_wt_spf_5xfad_comb) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))


# Extract feature of interest and check difference between strains colonized
exclusive_5xfad_liver <-  rbind(combined_unique_5xfad_5xfad_spf, combined_unique_5xfad_wt_spf)

data_5xfad_feat <- data_5xfad %>% dplyr::filter(!(SampleID %in% samples_to_remove$SampleID)) %>%
  dplyr::filter(!(SampleID %in% c("L5xfad_134", "L5xfad_98"))) %>%
  dplyr::select("SampleID", combined_unique_5xfad_5xfad_spf$ID, combined_unique_5xfad_wt_spf$ID) %>%
  dplyr::mutate(WT = rowSums(select(., combined_unique_5xfad_wt_spf$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., combined_unique_5xfad_5xfad_spf$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_metabolomics)

vip_5xfad_spf_wt <- VIPs_5xfad_spf_strain_Load %>% dplyr::filter(GroupContrib == "WT")
vip_5xfad_spf_mut <- VIPs_5xfad_spf_strain_Load %>% dplyr::filter(GroupContrib == "Mut")

data_5xfad_vip_spf <- data_5xfad %>% dplyr::filter(!(SampleID %in% samples_to_remove$SampleID)) %>%
  dplyr::filter(!(SampleID %in% c("L5xfad_134", "L5xfad_98"))) %>%
  dplyr::select("SampleID", VIPs_5xfad_spf_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_5xfad_spf_wt$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_5xfad_spf_mut$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_metabolomics)

# Plot ratio
data_5xfad_feat %>% 
  dplyr::filter(Colonized == "SPF") %>%
  ggboxplot(y = "Ratio", x = "Strain", add = "jitter") + 
  stat_compare_means()

plot_ratio_liver_5xfad <- data_5xfad_vip_spf %>% 
  dplyr::filter(Colonized == "SPF") %>% 
  dplyr::mutate(Strain = factor(Strain, levels = c("WT", "Mut"))) %>%
  ggboxplot(y = "Ratio", x = "Strain", add = "jitter", ylab = "Ln(WT/Mut)",
            add.params = list(color = "Strain", alpha = 0.6), legend = "none",
            palette = c("#E69A8DFF", "#5F4B8BFF"), xlab = "Liver",
            title = "Differential features from PLS-DA models") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = plot_ratio_liver_5xfad, filename = "5xfad_liver_ratio.svg", device = "svg", dpi = "retina", width = 1.5, height = 2.5)
#write_csv(x = data_5xfad_vip_spf, file = "liver_ratio_5xfad.csv")

# Cliff's Delta
cliff.delta((data_5xfad_feat %>%
               dplyr::filter(Colonized == "SPF" & Strain == "WT"))$Ratio, 
            (data_5xfad_feat %>%
               dplyr::filter(Colonized == "SPF" & Strain == "Mut"))$Ratio)


# Extract for MN
liver_5xfad_interest <- VIPs_5xfad_spf_strain_Load %>% 
  dplyr::select(ID, comp1, GroupContrib) %>% 
  left_join(VIPs_5xfad_gf_strain_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>% 
  left_join(VIPs_5xfad_mut_colo_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>% 
  left_join(VIPs_5xfad_wt_colo_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>% 
  dplyr::mutate(VIP_mean = rowMeans(dplyr::select(., comp1.x, comp1.x.x, comp1.y, comp1.y.y), na.rm = TRUE)) %>%
  left_join(exclusive_5xfad_liver %>% dplyr::select(ID, Group)) %>% arrange(desc(VIP_mean))

colnames(liver_5xfad_interest) <- c("ID", "VIP1", "SPF_genotype", "VIP2", "GF_genotype",
                                    "VIP3", "5xFAD_colonization", "VIP4", "WT_colonization", 
                                    "VIP_mean", "Exclusivity")

liver_5xfad_interest <- liver_5xfad_interest %>%   
  dplyr::mutate(Genotype = case_when(SPF_genotype == GF_genotype ~ "Yes",
                                     TRUE ~ "No")) %>%
  dplyr::mutate(Microbiome = case_when(`5xFAD_colonization` == WT_colonization ~ "Yes",
                                       TRUE ~ "No")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))

#write_csv(x = liver_5xfad_interest, file = "5xfad_liver_features_info_mn.csv")


###########
# Check overlapping between 3xTG and 5xFAD
###########
overlap_liver <- liver_3xtg_interest %>% 
  inner_join(liver_5xfad_interest, by = "ID") %>%
  dplyr::mutate(Study_concordance = case_when(SPF_genotype.x == SPF_genotype.y ~ "Yes",
                                              TRUE ~ "No"))


#############
# Export table for Joint-RPCA
#############
data_3xtg_join <- data_3xtg %>%
  dplyr::left_join(metadata_metabolomics %>% dplyr::select(SampleID, sample_name)) %>%
  dplyr::relocate(sample_name, .after=SampleID) %>% 
  dplyr::mutate(sample_name = gsub(".Liver", "", sample_name)) %>% dplyr::select(-SampleID)

#write_csv(x = data_3xtg_join, file = "liver_3xtg_metabolomics_join.csv")

data_5xfad_join <- data_5xfad_filter %>%
  dplyr::left_join(metadata_metabolomics %>% dplyr::select(SampleID, sample_name)) %>%
  dplyr::relocate(sample_name, .after=SampleID) %>% 
  dplyr::mutate(sample_name = gsub(".Liver", "", sample_name)) %>% dplyr::select(-SampleID)

#write_csv(x = data_5xfad_join, file = "liver_5xfad_metabolomics_join.csv")


######
# Explore features that are commonly separating groups between two different studies
######
features_strains <- VIPs_3xtg_strain_Load %>% 
  inner_join(VIPs_5xfad_strain_Load %>% dplyr::select(ID, comp1, GroupContrib), by = "ID") %>%
  dplyr::mutate(Concordance = case_when(`GroupContrib.x` == `GroupContrib.y` ~ paste("Yes", `GroupContrib.y`, sep = "_"),
                                        TRUE ~ "No")) %>% arrange(desc(Concordance)) %>%
  dplyr::mutate(Sum_VIP = `comp1.x` + `comp1.y`)

features_colonization <- VIPs_3xtg_colo_Load %>% 
  inner_join(VIPs_5xfad_colo_Load %>% dplyr::select(ID, comp1, GroupContrib), by = "ID")  %>%
  dplyr::mutate(Concordance = case_when(`GroupContrib.x` == `GroupContrib.y` ~ paste("Yes", `GroupContrib.y`, sep = "_"),
                                        TRUE ~ "No")) %>% arrange(desc(Concordance)) %>%
  dplyr::mutate(Sum_VIP = `comp1.x` + `comp1.y`)

features_sex <- VIPs_3xtg_sex_Load %>% 
  inner_join(VIPs_5xfad_sex_Load %>% dplyr::select(ID, comp1, GroupContrib), by = "ID")  %>%
  dplyr::mutate(Concordance = case_when(`GroupContrib.x` == `GroupContrib.y` ~ paste("Yes", `GroupContrib.y`, sep = "_"),
                                        TRUE ~ "No")) %>% arrange(desc(Concordance)) %>%
  dplyr::mutate(Sum_VIP = `comp1.x` + `comp1.y`)

# Check exclusive features
features_strains_3xtg <- VIPs_3xtg_strain_Load %>% 
  anti_join(VIPs_5xfad_strain_Load %>% dplyr::select(ID, comp1, GroupContrib), by = "ID")

features_colonization_3xtg <- VIPs_3xtg_colo_Load %>%  
  anti_join(VIPs_5xfad_colo_Load %>% dplyr::select(ID, comp1, GroupContrib), by = "ID")

features_sex_3xtg <- VIPs_3xtg_sex_Load %>% 
  anti_join(VIPs_5xfad_sex_Load %>% dplyr::select(ID, comp1, GroupContrib), by = "ID") 

features_strains_5xfad <- VIPs_5xfad_strain_Load %>% 
  anti_join(VIPs_3xtg_strain_Load %>% dplyr::select(ID, comp1, GroupContrib), by = "ID")

features_colonization_5xfad <- VIPs_5xfad_colo_Load %>%  
  anti_join(VIPs_3xtg_colo_Load %>% dplyr::select(ID, comp1, GroupContrib), by = "ID")

features_sex_5xfad <- VIPs_5xfad_sex_Load %>% 
  anti_join(VIPs_3xtg_sex_Load %>% dplyr::select(ID, comp1, GroupContrib), by = "ID") 


##############
# Filter MGF #
##############

VIPs_combined <- rbind(VIPs_3xtg_strain_Load, VIPs_3xtg_colo_Load, VIPs_3xtg_sex_Load,
                       VIPs_5xfad_strain_Load, VIPs_5xfad_colo_Load, VIPs_5xfad_sex_Load) %>%
  distinct(ID, .keep_all = TRUE) %>% arrange(ID)

# Import full mgf and filter it
dda <- Spectra("data/mice_caltech/liver/gnps_liver.mgf", source = MsBackendMgf())
dda_ids <- data.frame(ID = dda@backend@spectraData@listData$FEATURE_ID) %>%
  dplyr::mutate(Interest = ID %in% VIPs_combined$ID)
dda_filtered <- dda[dda_ids$Interest]
#export(dda_filtered, MsBackendMgf(), file = "Liver_VIPs_microbemasst.mgf", exportTitle = FALSE)
