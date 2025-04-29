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
feature_table <- read_csv("data/mice_caltech/fecal_gf/gnps_quant_fecal_gf.csv")
metadata <- read_csv("data/mice_caltech/fecal_gf/meta_fecal_gf.csv")
sample_order <- read.csv("data/mice_caltech/fecal_gf/sequence_fecal_gf.csv")
annotations <- read.delim("data/mice_caltech/fecal_gf/fbmn/nf_output/library/merged_results_with_gnps.tsv") %>%
  dplyr::filter(!str_detect(pattern = "REFRAME", LibraryName)) # remove REFRAME library
annotations$X.Scan. <- as.character(annotations$X.Scan.)
canopus <- read_tsv("data/mice_caltech/fecal_gf/canopus_fecal_gf.tsv") %>% dplyr::select(1,3,5,7,9)
canopus$id <- gsub("^.*sirius_", "", canopus$id)

sample_order$Plate <- gsub("^(.*?):.*", "\\1", sample_order$Position)

info_feature <- feature_table %>% dplyr::select(1:3,7)
colnames(info_feature) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature$Feature <- as.character(info_feature$Feature)

info_feature_complete <- info_feature %>% 
  left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  dplyr::select(1:5,18,24) %>% 
  left_join(canopus %>% distinct(id, .keep_all = TRUE), by = c("Feature" = "id"))

#write_csv(x = info_feature_complete, file = "fecal_gf_annotations.csv")

# Extra bile acids libraries and rev cosine
ba_ipsita <- read_csv("data/mice_caltech/fecal_gf/merged_Bile_acid_classic_networking.csv") %>% 
  dplyr::filter(SpectrumID %in% annotations$SpectrumID)
ba_rev <- read_tsv("data/mice_caltech/fecal_gf/candidate_BA_with_best_annotations.tsv") %>% 
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
  dplyr::select(-c("sample_type")) %>%
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

# Some samples have been re-run. Remove bad ones

# Sample rerun
sample_rerun <- data %>% dplyr::filter(str_detect(pattern = "_2024", SampleID)) %>%
  dplyr::select(SampleID)

sample_rerun$SampleID1 <- sub("_[^_]*$", "", sample_rerun$SampleID)

data_filter <- data %>% dplyr::filter(!(SampleID %in% sample_rerun$SampleID1)) %>%
  dplyr::filter(!(SampleID %in% c("Pool_qc_feces_37", "Pool_qc_feces_38", "Pool_qc_feces_39", 
                                  "Pool_qc_feces_40", "Pool_qc_feces_41", "Pool_qc_feces_42")))
data_filter$SampleID <- sub("_2024.*", "", data_filter$SampleID)

# Investigate total peak area
data_TIC <- data.frame(TIC = rowSums(data_filter %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("File.Name") %>% left_join(sample_order)

data_TIC %>%
  dplyr::mutate(SampleType = case_when(str_detect(pattern = "3xtg|5xfad", File.Name) ~ "Sample",
                                                  str_detect(pattern = "Pool", File.Name) ~ "Pool",
                                                             str_detect(pattern = "sixmix", File.Name) ~ "SixMix",
                                                                        TRUE ~ "Other")) %>%
  ggscatter("Order", "TIC", add = "reg.line", color = "SampleType", legend = "right") +
  stat_cor()

# looks like that around injection 230 there is clear increase in TIC. Tank change?

# Check sample type
sample_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "3xtg|5xfad", File.Name)) %>% summarise(mean(TIC))
pool_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Pool", File.Name)) %>% summarise(mean(TIC))
six_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "sixmix", File.Name)) %>% summarise(mean(TIC))
blank_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", File.Name)) %>% summarise(mean(TIC))
SRM_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "srm", File.Name)) %>% summarise(mean(TIC))

# Check TIC overtime for QCpool and Blanks
data_TIC %>% dplyr::filter(str_detect(pattern = "Pool", File.Name)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 1.5e9) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "sixmix", File.Name)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 5e8) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", File.Name)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 5e8) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "srm", File.Name)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 4e8) +
  stat_cor()

# Check all QCs --> look like there is a clear difference in two part of the run
data_TIC_six <- data_TIC  %>% dplyr::filter(str_detect(pattern = "sixmix", File.Name)) %>% arrange(Order)
mean(data_TIC_six$TIC[19:39] / data_TIC_six$TIC[1:18]) #1.35

data_TIC_pool <- data_TIC  %>% dplyr::filter(str_detect(pattern = "Pool", File.Name)) %>% arrange(Order)
mean(data_TIC_six$TIC[21:35] / data_TIC_six$TIC[1:20]) #1.34

data_TIC_srm <- data_TIC  %>% dplyr::filter(str_detect(pattern = "srm", File.Name)) %>% arrange(Order)
mean(data_TIC_srm$TIC[7:14] / data_TIC_srm$TIC[1:6]) #1.42

# Ratio are very similar --> Check the internal standard in the samples

# Check internal standard
fbmn_IS <- annotations %>% dplyr::filter(str_detect(Compound_Name, regex("sulf", ignore_case = TRUE))) %>% 
  distinct(X.Scan., .keep_all = TRUE) %>% dplyr::filter(Organism != "BILELIB19")

# Extract IS features for the table
table_IS <- data_filter %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% 
  dplyr::filter(ID %in% fbmn_IS$X.Scan.) %>% column_to_rownames("ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SampleID") %>% dplyr::filter(!(str_detect(SampleID, "srm|sixmix|Blank|Pool|pool"))) %>%
  dplyr::select(SampleID, `8173`) %>% left_join(metadata_metabolomics)

colnames(table_IS)[2] <- "Sulfadimethoxine"

table_IS %>% ggscatter(x = "Order", y = "Sulfadimethoxine", add = c("reg.line")) + ylim(0, 6e6) +
  stat_cor()

# Clear separation between first and second part of the run!
# There are some outlier samples with low or high IS values. Will be removed

data_IS_sample <- table_IS %>% arrange(Order)
mean(data_IS_sample$Sulfadimethoxine[180:376] / data_IS_sample$Sulfadimethoxine[1:179]) 

# Again, ratio is very similar --> 1.38

# Multiply values first part of the run by obtained ratio to account TIC shift in the middle of the run
table_is_fix <- data_IS_sample %>%
  dplyr::mutate(Sulfadimethoxine_fix = case_when(Order %in% c(10:227) ~ Sulfadimethoxine*1.3, 
                                                 TRUE ~ Sulfadimethoxine))

# Remove outliers
table_is_fix_out <- table_is_fix %>% 
  dplyr::filter(Sulfadimethoxine_fix > 1500000) %>%
  dplyr::filter(Sulfadimethoxine_fix < 4000000) 

table_is_fix_out %>%
  ggscatter(x = "Order", y = "Sulfadimethoxine_fix", add = c("reg.line")) +
  stat_cor()

table_is_fix_out %>% ggbarplot(x = "Order", y = "Sulfadimethoxine_fix", xlab = "Run Order", 
                       ylab = "Peak Area Sulfadimethoxine", title = "Internal Standard Acquisition") +
  geom_hline(yintercept = mean(table_is_fix_out$Sulfadimethoxine_fix, na.rm = TRUE), linetype = "dashed", color = "blue")

cv_is <- sd(table_is_fix_out$Sulfadimethoxine_fix)/mean(table_is_fix_out$Sulfadimethoxine_fix)


# Check Sample, QCPool, QCmix, and SRM
data_TIC %>% dplyr::filter(str_detect(pattern = "3xtg|5xfad", File.Name)) %>%
  dplyr::mutate(TIC_fix = case_when(Order %in% c(10:227) ~ TIC*1.3, 
                                                 TRUE ~ TIC)) %>%
  ggscatter("Order", "TIC_fix", add = "reg.line") + ylim(0, 2e9) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "Pool", File.Name)) %>%
  dplyr::mutate(TIC_fix = case_when(Order %in% c(10:227) ~ TIC*1.3, 
                                    TRUE ~ TIC)) %>%
  ggscatter("Order", "TIC_fix", add = "reg.line") + ylim(0, 1.5e9) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "sixmix", File.Name)) %>%
  dplyr::mutate(TIC_fix = case_when(Order %in% c(10:227) ~ TIC*1.4, 
                                    TRUE ~ TIC)) %>%
  ggscatter("Order", "TIC_fix", add = "reg.line") + ylim(0, 5e8) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "srm", File.Name)) %>%
  dplyr::mutate(TIC_fix = case_when(Order %in% c(10:227) ~ TIC*1.4, 
                                    TRUE ~ TIC)) %>%
  ggscatter("Order", "TIC_fix", add = "reg.line") + ylim(0, 5e8) +
  stat_cor()


# Fix the whole dataset (first 230 samples)
data_to_fix <- data_filter %>% left_join(metadata_metabolomics %>% dplyr::select(SampleID, Order)) %>%
  arrange(Order) %>% dplyr::mutate_at("Order", as.numeric) %>% dplyr::filter(Order %in% c(10:227)) %>%
  column_to_rownames("SampleID") %>% dplyr::select(-Order) %>% mutate_all(~ . * 1.3) %>% rownames_to_column("SampleID")

data_not_to_fix <- data_filter %>% left_join(metadata_metabolomics %>% dplyr::select(SampleID, Order)) %>%
  arrange(Order) %>% dplyr::mutate_at("Order", as.numeric) %>% dplyr::filter(!(Order %in% c(10:227))) %>% 
  dplyr::select(-Order)

data_fix_final <- rbind(data_to_fix, data_not_to_fix)


# Check features per sample type
data_blank <- data_fix_final %>% dplyr::filter(str_detect(pattern = "Blank", SampleID))
data_pool <- data_fix_final %>% dplyr::filter(str_detect(pattern = "Pool", SampleID))
data_sixmix <- data_fix_final %>% dplyr::filter(str_detect(pattern = "sixmix", SampleID))
data_srm <- data_fix_final %>% dplyr::filter(str_detect(pattern = "srm", SampleID))

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
# CVs of the six standards are ranging between 2% to 8% --> very good

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
data_clean <- data_fix_final %>% dplyr::select(-c(feature_to_remove$Feature)) %>%
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
feature_to_remove_rt <- info_feature_complete %>% dplyr::filter(RT < 0.2 | RT > 8) %>%
  dplyr::filter(!(Feature %in% feature_to_remove$Feature)) %>%
  dplyr::filter(!(Feature %in% feature_to_remove_mix$Feature))

# Final cleaned table
data_clean3 <- data_clean2 %>% dplyr::select(-c(feature_to_remove_rt$Feature))

# PCA raw data
PCA_raw <- mixOmics::pca(data_clean3 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

i <- "intervention"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Fecal", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Keep only samples
data_sample <- data_clean3 %>% 
  dplyr::filter(!(str_detect(pattern = "srm|sixmix|Pool|Blank|blank|pool", SampleID))) %>% 
  dplyr::filter(!(SampleID %in% c("feces_3xtg_70", "feces_5xfad_61", "feces_5xfad_62")))

# PCA sample data
PCA_raw <- mixOmics::pca(data_sample %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

i <- "intervention"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Fecal ", i, sep = " "),
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

i <- "intervention"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = "Fecal - GF", palette = c("#B7E6A5", "#007188"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = PCA_plot, filename = "fecal_gf_pca_metabolome.svg", device = "svg", dpi = "retina", height = 2, width = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_sample_clr, method = "euclidean")
disper_intervention <- betadisper(dist_metabolites, PCA_whole_scores$intervention)
anova(disper_intervention)
permanova <- adonis2(dist_metabolites ~ intervention, PCA_whole_scores, na.action = na.omit)

# Clear separation between sacrifice and longitudinal data

# Separate animal by design, study, colonization and strain
sample_long <- metadata_metabolomics %>% dplyr::filter(intervention == "Longitudinal")
sample_sac <- metadata_metabolomics %>% dplyr::filter(intervention == "Sacrifice")

sample_3xtg <- metadata_metabolomics %>% dplyr::filter(Study == "3xTG")
sample_5xfad <- metadata_metabolomics %>% dplyr::filter(Study == "5xFAD")

sample_wt <- metadata_metabolomics %>% dplyr::filter(Strain == "WT")
sample_mut <- metadata_metabolomics %>% dplyr::filter(Strain == "Mut")

sample_male <- metadata_metabolomics %>% dplyr::filter(sex == "male")
sample_female <- metadata_metabolomics %>% dplyr::filter(sex == "female")


##################
# SACRIFICE DATA #
##################
data_sac <- data_sample %>% dplyr::filter(SampleID %in% sample_sac$SampleID)

# RCLR transformation
data_sacrifice_clr <- decostand(data_sac %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sacrifice_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PCA_sacrifice_plots <- list()

for (i in c("box", "host_age", "sex", "mouse_strain", "Plate", "Strain", "Study")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Fecal RCLR Sacrifice", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_sacrifice_plots[[i]] <- PCA_plot
  
}


PCA_sacrifice_plots_final <- wrap_plots(PCA_sacrifice_plots, nrow = 3)

# PERMANOVA
dist_metabolites <- vegdist(data_sacrifice_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$mouse_strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ mouse_strain * sex + 
                       box + Plate, PCA_whole_scores, na.action = na.omit)

i <- "Study"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = "Fecal GF - Sacrifice", palette = c("#009392", "#F0746E"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = PCA_plot, filename = "fecal_gf_sacrifice_pca_metabolome.svg", device = "svg", dpi = "retina", height = 2, width = 2)

#write_csv(x = data_sac, file = "fecal_gf_sacrifice_metabolome.csv")
#write_csv(x = data_sacrifice_clr, file = "fecal_gf_sacrifice_metabolome_rclr.csv")
#write_csv(x = metadata_metabolomics %>% dplyr::filter(SampleID %in% data_sac$SampleID), file = "fecal_gf_sacrifice_metadata.csv")


##################
# Sacrifice 3xTG #
##################
data_sac_3xtg <- data_sac %>% dplyr::filter(SampleID %in% sample_3xtg$SampleID)

# RCLR transformation
data_sac_3xtg_clr <- decostand(data_sac_3xtg %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sac_3xtg_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PCA_3xtg_sacrifice_plots <- list()

for (i in c("box", "host_age", "sex", "Strain")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Fecal RCLR 3xTG Sacrifice", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_3xtg_sacrifice_plots[[i]] <- PCA_plot
  
}

PCA_3xtg_sacrifice_plots_final <- wrap_plots(PCA_3xtg_sacrifice_plots, nrow = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_sac_3xtg_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$Strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ Strain * sex + 
                       box + host_age, PCA_whole_scores, na.action = na.omit)


# PLSDA Sacrificed 3xTG - Strain
PLSDA_3xtg_sac_strain <- mixOmics::plsda(data_sac_3xtg_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                         PCA_whole_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_3xtg_sac_strain_scores <- data.frame(PLSDA_3xtg_sac_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_strain_plot <- PLSDA_3xtg_sac_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Fecal 3xTG Sacrificed - Strain",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_sac_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_sac_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_sac_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_sac_strain <- plotLoadings(PLSDA_3xtg_sac_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_sac_strain <- perf(PLSDA_3xtg_sac_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_sac_strain, legend = FALSE)

VIPs_3xtg_sac_strain <- as.data.frame(mixOmics::vip(PLSDA_3xtg_sac_strain))
VIPs_3xtg_sac_strain_filter <- dplyr::filter(VIPs_3xtg_sac_strain, VIPs_3xtg_sac_strain$comp1 > 1)
VIPs_3xtg_sac_strain_filter$ID <- rownames(VIPs_3xtg_sac_strain_filter)
VIPs_3xtg_sac_strain_select <- VIPs_3xtg_sac_strain_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_sac_strain_Load <- VIPs_3xtg_sac_strain_select %>% 
  left_join(Loadings_3xtg_sac_strain, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))  %>% arrange(desc(comp1))

# PLSDA Sacrificed 3xTG - Sex
PLSDA_3xtg_sac_sex <- mixOmics::plsda(data_sac_3xtg_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                      PCA_whole_scores$sex, ncomp = 2, scale = TRUE)
PLSDA_3xtg_sac_sex_scores <- data.frame(PLSDA_3xtg_sac_sex$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_sac_sex_plot <- PLSDA_3xtg_sac_sex_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "sex", alpha = 0.6, title = "PLSDA - Fecal 3xTG Sacrificed - Sex",
            xlab = paste("Component 1 (", round(PLSDA_3xtg_sac_sex$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_3xtg_sac_sex$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_3xtg_sac_sex_scores %>% group_by(sex) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = sex), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_3xtg_sac_sex <- plotLoadings(PLSDA_3xtg_sac_sex, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_3xtg_sac_sex <- perf(PLSDA_3xtg_sac_sex, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_3xtg_sac_sex, legend = FALSE)

VIPs_3xtg_sac_sex <- as.data.frame(mixOmics::vip(PLSDA_3xtg_sac_sex))
VIPs_3xtg_sac_sex_filter <- dplyr::filter(VIPs_3xtg_sac_sex, VIPs_3xtg_sac_sex$comp1 > 1)
VIPs_3xtg_sac_sex_filter$ID <- rownames(VIPs_3xtg_sac_sex_filter)
VIPs_3xtg_sac_sex_select <- VIPs_3xtg_sac_sex_filter %>% dplyr::select(ID, comp1)
VIPs_3xtg_sac_sex_Load <- VIPs_3xtg_sac_sex_select %>% 
  left_join(Loadings_3xtg_sac_sex, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# Remove feature influenced by sex
VIPs_3xtg_sac_strain_no_sex <- VIPs_3xtg_sac_strain_Load %>% 
  anti_join(VIPs_3xtg_sac_sex_Load %>% dplyr::select(ID))

#write_csv(x = VIPs_3xtg_sac_strain_no_sex, "3xtg_fecal_gf_plsda_vip_all_samples_genotype.csv")


############
# Compare 3xtg male v wt male to 3xtg female v wt female to see genotype effect
############
data_3xtg_male <- data_sac_3xtg %>% 
  dplyr::filter(SampleID %in% sample_male$SampleID)

data_3xtg_female <- data_sac_3xtg %>% 
  dplyr::filter(SampleID %in% sample_female$SampleID)

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
              title = paste("PCA - Brain 3xTG Male", i, sep = " "),
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
              title = paste("PCA - Brain 3xTG Female", i, sep = " "),
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
                                          PCA_3xtg_male_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_3xtg_male_strain_scores <- data.frame(PLSDA_3xtg_male_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_male_strain_plot <- PLSDA_3xtg_male_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Fecal GF 3xTG Male - Strain",
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

# PLSDA 3xTG Female - Strain
PLSDA_3xtg_female_strain <- mixOmics::plsda(data_3xtg_female_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                            PCA_3xtg_female_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_3xtg_female_strain_scores <- data.frame(PLSDA_3xtg_female_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_3xtg_female_strain_plot <- PLSDA_3xtg_female_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Fecal GF 3xTG Female - Strain",
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

# Sample matching 
list_3xtg_sex <- list(
  `Male 3xTG` = (VIPs_3xtg_male_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `Male WT` =  (VIPs_3xtg_male_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID,
  `Female 3xTG` = (VIPs_3xtg_female_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `Female WT` = (VIPs_3xtg_female_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID)

sex_3xtg <- UpSetR::upset(fromList(list_3xtg_sex), nsets = 4, nintersects = 6, 
                             point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                             sets = c("Male 3xTG", "Male WT", "Female 3xTG", "Female WT"),
                             queries = list(list(query = intersects, params = list("Female 3xTG", 
                                                                                   "Male 3xTG"), color = "steelblue1", active = T), 
                                            list(query = intersects, params = list("Female WT", 
                                                                                   "Male WT"),color = "orange", active = T)))

#pdf(file="3xtg_fecal_gf_upset_genotype.pdf", width = 3, height = 2.5)
#sex_3xtg
#dev.off()


# Check features of interest
vip_3xtg_gf_wt <- VIPs_3xtg_sac_strain_Load %>% dplyr::filter(GroupContrib == "WT")
vip_3xtg_gf_mut <- VIPs_3xtg_sac_strain_Load %>% dplyr::filter(GroupContrib == "Mut")

data_3xtg_vip_spf <- data_sac_3xtg %>%
  dplyr::select("SampleID", VIPs_3xtg_sac_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_3xtg_gf_wt$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_3xtg_gf_mut$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_metabolomics)

plot_ratio <- data_3xtg_vip_spf %>%
  dplyr::mutate(Strain = factor(Strain, levels = c("WT", "Mut"))) %>%
  ggboxplot(y = "Ratio", x = "Strain", add = "jitter", title = "Fecal GF - 3xTG", 
            ylab = "Log(Top WT/Top MUT)", xlab = FALSE) + 
  stat_compare_means()+
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8)) # color by sex and age to see an effect

#ggsave(plot = plot_ratio, filename = "3xtg_fecal_gf_ratio.svg", device = "svg", dpi = "retina", width = 2, height = 3)

# Cliff's Delta
cliff.delta((data_3xtg_vip_spf %>%
               dplyr::filter(Strain == "WT"))$Ratio, 
            (data_3xtg_vip_spf %>%
               dplyr::filter(Strain == "Mut"))$Ratio)


# Extract features of interest for MN
fecal_gf_3xtg_interest <- VIPs_3xtg_sac_strain_Load %>% 
  dplyr::select(ID, comp1, GroupContrib) %>% 
  left_join(VIPs_3xtg_male_strain_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>% 
  left_join(VIPs_3xtg_female_strain_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>%
  dplyr::mutate(VIP_mean = rowMeans(dplyr::select(., comp1.x, comp1.y, comp1), na.rm = TRUE))

colnames(fecal_gf_3xtg_interest) <- c("ID", "VIP1", "All_genotype", "VIP2", "Male_genotype",
                                   "VIP3", "Female_genotype", "VIP_mean")

fecal_gf_3xtg_interest <- fecal_gf_3xtg_interest %>%   
  dplyr::mutate(Genotype = case_when(Male_genotype == Female_genotype ~ "Yes",
                                     TRUE ~ "No")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))

#write_csv(x = fecal_gf_3xtg_interest, file = "3xtg_fecal_gf_features_info_mn.csv")


###################
# Sacrifice 5xFAD #
###################
data_sac_5xfad <- data_sac %>% 
  dplyr::filter(SampleID %in% sample_5xfad$SampleID)

# RCLR transformation
data_sac_5xfad_clr <- decostand(data_sac_5xfad %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sac_5xfad_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PCA_5xfad_sacrifice_plots <- list()

for (i in c("box", "host_age", "sex", "Strain")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Fecal RCLR 5xFAD Sacrifice", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_5xfad_sacrifice_plots[[i]] <- PCA_plot
  
}

PCA_5xfad_sacrifice_plots_final <- wrap_plots(PCA_5xfad_sacrifice_plots, nrow = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_sac_5xfad_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$mouse_strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ Strain * sex + box + 
                       Plate, PCA_whole_scores, na.action = na.omit)

# looks like a bunch of samples is unusable. remove them
sample_sac_5xfad_to_remove <- PCA_whole_scores %>% dplyr::filter(PC1 > 30)

data_sac_5xfad_filter <- data_sac_5xfad %>% 
  dplyr::filter(!(SampleID %in% sample_sac_5xfad_to_remove$SampleID))

# RCLR transformation
data_sac_5xfad_clr <- decostand(data_sac_5xfad_filter %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sac_5xfad_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PCA_5xfad_sacrifice_plots <- list()

for (i in c("box", "host_age", "sex", "Strain")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Fecal RCLR 5xFAD Sacrifice", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_5xfad_sacrifice_plots[[i]] <- PCA_plot
  
}

PCA_5xfad_sacrifice_plots_final <- wrap_plots(PCA_5xfad_sacrifice_plots, nrow = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_sac_5xfad_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$mouse_strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ Strain * sex + box + 
                       Plate, PCA_whole_scores, na.action = na.omit)


# PLSDA Sacrificed 5xFAD - Strain
PLSDA_5xfad_sac_strain <- mixOmics::plsda(data_sac_5xfad_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                          PCA_whole_scores$Strain, ncomp = 3, scale = TRUE)
PLSDA_5xfad_sac_strain_scores <- data.frame(PLSDA_5xfad_sac_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_sac_strain_plot <- PLSDA_5xfad_sac_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Fecal 5xFAD Sacrificed - Strain",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_sac_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_sac_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_sac_strain_scores %>% group_by(Strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_sac_strain <- plotLoadings(PLSDA_5xfad_sac_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_sac_strain <- perf(PLSDA_5xfad_sac_strain, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_sac_strain, legend = FALSE)

VIPs_5xfad_sac_strain <- as.data.frame(mixOmics::vip(PLSDA_5xfad_sac_strain))
VIPs_5xfad_sac_strain_filter <- dplyr::filter(VIPs_5xfad_sac_strain, VIPs_5xfad_sac_strain$comp1 > 1)
VIPs_5xfad_sac_strain_filter$ID <- rownames(VIPs_5xfad_sac_strain_filter)
VIPs_5xfad_sac_strain_select <- VIPs_5xfad_sac_strain_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_sac_strain_Load <- VIPs_5xfad_sac_strain_select %>% 
  left_join(Loadings_5xfad_sac_strain, by = c("ID" = "rowname")) %>% arrange(desc(comp1)) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))

# PLSDA Sacrificed 5xFAD - Sex
PLSDA_5xfad_sac_sex <- mixOmics::plsda(data_sac_5xfad_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                       PCA_whole_scores$sex, ncomp = 2, scale = TRUE)
PLSDA_5xfad_sac_sex_scores <- data.frame(PLSDA_5xfad_sac_sex$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_sac_sex_plot <- PLSDA_5xfad_sac_sex_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "sex", alpha = 0.6, title = "PLSDA - Fecal 5xFAD Sacrificed - Sex",
            xlab = paste("Component 1 (", round(PLSDA_5xfad_sac_sex$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_5xfad_sac_sex$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_5xfad_sac_sex_scores %>% group_by(sex) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = sex), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_5xfad_sac_sex <- plotLoadings(PLSDA_5xfad_sac_sex, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_5xfad_sac_sex <- perf(PLSDA_5xfad_sac_sex, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_5xfad_sac_sex, legend = FALSE)

VIPs_5xfad_sac_sex <- as.data.frame(mixOmics::vip(PLSDA_5xfad_sac_sex))
VIPs_5xfad_sac_sex_filter <- dplyr::filter(VIPs_5xfad_sac_sex, VIPs_5xfad_sac_sex$comp1 > 1)
VIPs_5xfad_sac_sex_filter$ID <- rownames(VIPs_5xfad_sac_sex_filter)
VIPs_5xfad_sac_sex_select <- VIPs_5xfad_sac_sex_filter %>% dplyr::select(ID, comp1)
VIPs_5xfad_sac_sex_Load <- VIPs_5xfad_sac_sex_select %>% 
  left_join(Loadings_5xfad_sac_sex, by = c("ID" = "rowname")) %>% arrange(desc(comp1)) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))


# Remove feature influenced by sex
VIPs_5xfad_sac_strain_no_sex <- VIPs_5xfad_sac_strain_Load %>% 
  anti_join(VIPs_5xfad_sac_sex_Load %>% dplyr::select(ID))

#write_csv(x = VIPs_5xfad_sac_strain_no_sex, file = "5xfad_fecal_gf_plsda_vip_all_samples_genotype.csv")


############
# Compare 5xfad male v wt male to 5xfad female v wt female to see strain and microbiome effect within sex
############
data_5xfad_male <- data_sac_5xfad_filter %>% 
  dplyr::filter(SampleID %in% sample_male$SampleID)

data_5xfad_female <- data_sac_5xfad_filter %>% 
  dplyr::filter(SampleID %in% sample_female$SampleID)

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
              title = paste("PCA - Brain 5xFAD Male", i, sep = " "),
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
              title = paste("PCA - Brain 5xFAD Female", i, sep = " "),
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
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Fecal GF 5xFAD Male - Strain",
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

# PLSDA 5xFAD Female - Strain
PLSDA_5xfad_female_strain <- mixOmics::plsda(data_5xfad_female_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                             PCA_5xfad_female_scores$Strain, ncomp = 2, scale = TRUE)
PLSDA_5xfad_female_strain_scores <- data.frame(PLSDA_5xfad_female_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_5xfad_female_strain_plot <- PLSDA_5xfad_female_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Strain", alpha = 0.6, title = "PLSDA - Fecal GF 5xFAD Female - Strain",
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


# Sample matching 
list_5xfad_sex <- list(
  `Male 5xFAD` = (VIPs_5xfad_male_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `Male WT` =  (VIPs_5xfad_male_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID,
  `Female 5xFAD` = (VIPs_5xfad_female_strain_Load %>% dplyr::filter(GroupContrib == "Mut"))$ID,
  `Female WT` = (VIPs_5xfad_female_strain_Load %>% dplyr::filter(GroupContrib == "WT"))$ID)

sex_5xfad <- UpSetR::upset(fromList(list_5xfad_sex), nsets = 4, nintersects = 7, 
                          point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                          sets = c("Male 5xFAD", "Male WT", "Female 5xFAD", "Female WT"),
                          queries = list(list(query = intersects, params = list("Female 5xFAD", 
                                                                                "Male 5xFAD"), color = "steelblue1", active = T), 
                                         list(query = intersects, params = list("Female WT", 
                                                                                "Male WT"),color = "orange", active = T)))

#pdf(file="5xfad_fecal_gf_upset_genotype.pdf", width = 3, height = 2.5)
#sex_5xfad
#dev.off()

# Check features of interest
vip_5xfad_gf_wt <- VIPs_5xfad_sac_strain_Load %>% dplyr::filter(GroupContrib == "WT")
vip_5xfad_gf_mut <- VIPs_5xfad_sac_strain_Load %>% dplyr::filter(GroupContrib == "Mut")

data_5xfad_vip_spf <- data_sac_5xfad_filter %>%
  dplyr::select("SampleID", VIPs_5xfad_sac_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_5xfad_gf_wt$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_5xfad_gf_mut$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_metabolomics)

plot_ratio <- data_5xfad_vip_spf %>%
  dplyr::mutate(Strain = factor(Strain, levels = c("WT", "Mut"))) %>%
  ggboxplot(y = "Ratio", x = "Strain", add = "jitter", title = "Fecal GF - 5xFAD", 
            ylab = "Log(Top WT/Top MUT)", xlab = FALSE) + 
  stat_compare_means()+
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8)) # color by sex and age to see an effect

#ggsave(plot = plot_ratio, filename = "5xfad_fecal_gf_ratio.svg", device = "svg", dpi = "retina", width = 2, height = 3)

# Cliff's Delta
cliff.delta((data_5xfad_vip_spf %>%
               dplyr::filter(Strain == "WT"))$Ratio, 
            (data_5xfad_vip_spf %>%
               dplyr::filter(Strain == "Mut"))$Ratio)


# Extract features of interest for MN
fecal_gf_5xfad_interest <- VIPs_5xfad_sac_strain_Load %>% 
  dplyr::select(ID, comp1, GroupContrib) %>% 
  left_join(VIPs_5xfad_male_strain_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>% 
  left_join(VIPs_5xfad_female_strain_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>%
  dplyr::mutate(VIP_mean = rowMeans(dplyr::select(., comp1.x, comp1.y, comp1), na.rm = TRUE))

colnames(fecal_gf_5xfad_interest) <- c("ID", "VIP1", "All_genotype", "VIP2", "Male_genotype",
                                      "VIP3", "Female_genotype", "VIP_mean")

fecal_gf_5xfad_interest <- fecal_gf_5xfad_interest %>%   
  dplyr::mutate(Genotype = case_when(Male_genotype == Female_genotype ~ "Yes",
                                     TRUE ~ "No")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))

#write_csv(x = fecal_gf_5xfad_interest, file = "5xfad_fecal_gf_features_info_mn.csv")


###########
# Check overlapping between 3xTG and 5xFAD
###########
overlap_fecal_gf <- fecal_gf_3xtg_interest %>% 
  inner_join(fecal_gf_5xfad_interest, by = "ID") %>%
  dplyr::mutate(Study_concordance = case_when(All_genotype.x == All_genotype.y ~ "Yes",
                                              TRUE ~ "No"))


#############
# Export table for Joint-RPCA
#############
data_3xtg_join <- data_sac_3xtg %>%
  dplyr::left_join(metadata_metabolomics %>% dplyr::select(SampleID, sample_name)) %>%
  dplyr::relocate(sample_name, .after=SampleID) %>% 
  dplyr::mutate(sample_name = gsub("\\.[^.]*$", "", sample_name)) %>% dplyr::select(-SampleID)

#write_csv(x = data_3xtg_join, file = "fecal_gf_sac_3xtg_metabolomics_join.csv")

data_5xfad_join <- data_sac_5xfad %>%
  dplyr::left_join(metadata_metabolomics %>% dplyr::select(SampleID, sample_name)) %>%
  dplyr::relocate(sample_name, .after=SampleID) %>% 
  dplyr::mutate(sample_name = gsub("\\.[^.]*$", "", sample_name)) %>% dplyr::select(-SampleID)

#write_csv(x = data_5xfad_join, file = "fecal_gf_sac_5xfad_metabolomics_join.csv")


#####################
# LONGITUDINAL DATA #
#####################
data_long <- data_sample %>% 
  dplyr::filter(SampleID %in% sample_long$SampleID)

# RCLR transformation
data_long_clr <- decostand(data_long %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_long_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PCA_long_plots <- list()

for (i in c("box", "host_age", "host_subject_id", "sex", "mouse_strain", "Study", "Strain")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Longitudinal", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_long_plots[[i]] <- PCA_plot
  
}

PCA_long_plots_final <- wrap_plots(PCA_long_plots, nrow = 3)

# PERMANOVA
dist_metabolites <- vegdist(data_long_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$mouse_strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ mouse_strain * sex * host_age + 
                       box, PCA_whole_scores, na.action = na.omit)


#####################
# LONGITUDINAL 3xTG #
#####################
data_long_3xtg <- data_long %>% 
  dplyr::filter(SampleID %in% sample_3xtg$SampleID)

metadata_long_3xtg <- metadata_metabolomics %>% 
  dplyr::filter(intervention == "Longitudinal") %>%
  dplyr::filter(Study == "3xTG") %>%
  dplyr::mutate_at("host_age", as.numeric) %>% arrange(host_age)

# Plot longitudinal metadata information
meta_long_3xtg_plot <- metadata_long_3xtg %>% 
  group_by(host_age, mouse_strain) %>%   
  summarise(count = n(), .groups = "drop") %>%
  ggplot(aes(x = factor(host_age), y = count, fill = mouse_strain)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#5F4B8BFF", "#E69A8DFF")) +
  labs(x = "Age (months)", y = "Count", fill = "Genotype",
       title = "Fecal samples in 3xTg GF") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = meta_long_3xtg_plot, filename = "Metadata_fecal_3xtg_long.svg", device = "svg", dpi = "retina", width = 2.5, height = 2)


# RCLR transformation
data_long_3xtg_clr <- decostand(data_long_3xtg %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_long_3xtg_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_long_3xtg_plots <- list()

for (i in c("box", "host_age", "host_subject_id", "sex", "mouse_strain")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Longitudinal", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_long_3xtg_plots[[i]] <- PCA_plot
  
}

PCA_long_3xtg_plots_final <- wrap_plots(PCA_long_3xtg_plots, nrow = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_long_3xtg_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$Strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ Strain + host_age + sex + host_subject_id, PCA_whole_scores, na.action = na.omit)

i <- "Strain"

PCA_plot_3xtg_long <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, legend = "none",
            title = paste("PCA - Longitudinal", i, sep = " "), palette = c("#5F4B8BFF", "#E69A8DFF"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = PCA_plot_3xtg_long, filename = "PCA_3xtg_long.svg", device = "svg", dpi = "retina", height = 2, width = 2)


# Extract features and generate log ratio which was established in the Sacrifice cohort
data_3xtg_vip_gf_long <- data_long_3xtg %>%
  dplyr::select("SampleID", VIPs_3xtg_sac_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_3xtg_gf_wt$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_3xtg_gf_mut$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_long_3xtg)

plot_ratio_3xtg_long <- data_3xtg_vip_gf_long %>% 
  dplyr::mutate_at("host_age", as.numeric) %>%
  ggscatter(x = "host_age", y = "Ratio", color = "mouse_strain", add = "loess", ylab = "Ln(WTTop/MutTop)",
            xlab = "Age (months)", legend = "none", palette = c("#5F4B8BFF", "#E69A8DFF"), alpha = 0.4,
            title = "Natural Log Ratio Across Time in 3xTg Study") + 
  scale_x_continuous(breaks = seq(min(data_3xtg_vip_gf_long$host_age, na.rm = TRUE),
                                  max(data_3xtg_vip_gf_long$host_age, na.rm = TRUE), by = 1)) +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = plot_ratio_3xtg_long, filename = "Log_ratio_3xtg_gf_long.svg", device = "svg", dpi = "retina", width = 3.5, height = 2)


# linear mixed effect model to test separation
library(lmerTest)
model <- lmer(Ratio ~ host_age * Strain + (1 | host_subject_id), data = data_3xtg_vip_gf_long)
summary(model)

# check bile acid ratio only
vip_3xtg_gf_mut_ba <- vip_3xtg_gf_mut %>% 
  dplyr::filter(str_detect(pattern = "butanoic acid|pentanoic acid|bile acid|dihydroxyoxane-2-carboxylic acid", Compound_Name))
vip_3xtg_gf_wt_ba <- vip_3xtg_gf_wt %>% 
  dplyr::filter(str_detect(pattern = "butanoic acid|pentanoic acid|bile acid|dihydroxyoxane-2-carboxylic acid", Compound_Name))

data_3xtg_vip_gf_long_ba <- data_long_3xtg %>%
  dplyr::select("SampleID", VIPs_3xtg_sac_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_3xtg_gf_wt_ba$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_3xtg_gf_mut_ba$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_long_3xtg)

plot_ratio_3xtg_long_ba <- data_3xtg_vip_gf_long_ba %>% 
  dplyr::mutate_at("host_age", as.numeric) %>%
  ggscatter(x = "host_age", y = "Ratio", color = "Strain", add = "loess", ylab = "Ln(WTTop/MutTop)",
            xlab = "Age (months)", legend = "none", palette = c("#5F4B8BFF", "#E69A8DFF"), alpha = 0.4,
            title = "Bile Acid Ratio in Fecal 3xTg GF") + 
  scale_x_continuous(breaks = seq(min(data_3xtg_vip_gf_long_ba$host_age, na.rm = TRUE),
                                  max(data_3xtg_vip_gf_long_ba$host_age, na.rm = TRUE), by = 1)) +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = plot_ratio_3xtg_long_ba, filename = "Log_ratio_3xtg_gf_long_ba.svg", device = "svg", dpi = "retina", width = 3.5, height = 2)

# linear mixed effect model to test separation
model_ba <- lmer(Ratio ~ host_age * Strain + (1 | host_subject_id), data = data_3xtg_vip_gf_long_ba)
summary(model_ba)


######################
# LONGITUDINAL 5xFAD #
######################
data_long_5xfad <- data_long %>% 
  dplyr::filter(SampleID %in% sample_5xfad$SampleID)

metadata_long_5xfad <- metadata_metabolomics %>% 
  dplyr::filter(intervention == "Longitudinal") %>%
  dplyr::filter(Study == "5xFAD") %>%
  dplyr::mutate_at("host_age", as.numeric) %>% arrange(host_age)

# Plot longitudinal metadata information
meta_long_5xfad_plot <- metadata_long_5xfad %>% 
  group_by(host_age, mouse_strain) %>%   
  summarise(count = n(), .groups = "drop") %>%
  ggplot(aes(x = factor(host_age), y = count, fill = mouse_strain)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#5F4B8BFF", "#E69A8DFF")) +
  labs(x = "Age (months)", y = "Count", fill = "Genotype",
       title = "Fecal samples in 5xFAD GF") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = meta_long_5xfad_plot, filename = "Metadata_fecal_5xfad_gf_long.svg", device = "svg", dpi = "retina", width = 2.5, height = 2)


# RCLR transformation
data_long_5xfad_clr <- decostand(data_long_5xfad %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_long_5xfad_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_long_5xfad_plots <- list()

for (i in c("box", "host_age", "host_subject_id", "sex", "mouse_strain")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Longitudinal", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_long_5xfad_plots[[i]] <- PCA_plot
  
}

PCA_long_5xfad_plots_final <- wrap_plots(PCA_long_5xfad_plots, nrow = 2)

i <- "Strain"

PCA_plot_5xfad_long <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, legend = "none",
            title = paste("PCA - Longitudinal", i, sep = " "), palette = c("#5F4B8BFF", "#E69A8DFF"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = PCA_plot_5xfad_long, filename = "PCA_5xfad_long.svg", device = "svg", dpi = "retina", height = 2, width = 2)


# PERMANOVA
dist_metabolites <- vegdist(data_long_5xfad_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$Strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ Strain + sex + host_age + host_subject_id, PCA_whole_scores, na.action = na.omit)


# Extract features and generate log ratio which was established in the Sacrifice cohort
data_5xfad_vip_gf_long <- data_long_5xfad %>%
  dplyr::select("SampleID", VIPs_5xfad_sac_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_5xfad_gf_wt$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_5xfad_gf_mut$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_long_5xfad)

plot_ratio_5xfad_long <- data_5xfad_vip_gf_long %>% 
  dplyr::mutate_at("host_age", as.numeric) %>%
  ggscatter(x = "host_age", y = "Ratio", color = "mouse_strain", add = "loess", ylab = "Ln(WTTop/MutTop)",
            xlab = "Age (months)", legend = "none", palette = c("#5F4B8BFF", "#E69A8DFF"), alpha = 0.4,
            title = "Natural Log Ratio Across Time in 5xFAD Study") + 
  scale_x_continuous(breaks = seq(min(data_5xfad_vip_gf_long$host_age, na.rm = TRUE),
                                  max(data_5xfad_vip_gf_long$host_age, na.rm = TRUE), by = 1)) +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = plot_ratio_5xfad_long, filename = "Log_ratio_5xfad_gf_long.svg",device = "svg", dpi = "retina", width = 3.5, height = 2)


# linear mixed effect model to test separation
library(lmerTest)
model <- lmer(Ratio ~ host_age * Strain + (1 | host_subject_id), data = data_5xfad_vip_gf_long)
summary(model)


# check bile acid ratio only
vip_5xfad_gf_mut_ba <- vip_5xfad_gf_mut %>% 
  dplyr::filter(str_detect(pattern = "butanoic acid|pentanoic acid|bile acid|dihydroxyoxane-2-carboxylic acid", Compound_Name))
vip_5xfad_gf_wt_ba <- vip_5xfad_gf_wt %>% 
  dplyr::filter(str_detect(pattern = "butanoic acid|pentanoic acid|bile acid|dihydroxyoxane-2-carboxylic acid", Compound_Name))

data_5xfad_vip_gf_long_ba <- data_long_5xfad %>%
  dplyr::select("SampleID", VIPs_5xfad_sac_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_5xfad_gf_wt_ba$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_5xfad_gf_mut_ba$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_long_5xfad)

plot_ratio_5xfad_long_ba <- data_5xfad_vip_gf_long_ba %>% 
  dplyr::mutate_at("host_age", as.numeric) %>%
  ggscatter(x = "host_age", y = "Ratio", color = "Strain", add = "loess", ylab = "Ln(WTTop/MutTop)",
            xlab = "Age (months)", legend = "none", palette = c("#5F4B8BFF", "#E69A8DFF"), alpha = 0.4,
            title = "Bile Acid Ratio in Fecal 5xfad GF") + 
  scale_x_continuous(breaks = seq(min(data_5xfad_vip_gf_long_ba$host_age, na.rm = TRUE),
                                  max(data_5xfad_vip_gf_long_ba$host_age, na.rm = TRUE), by = 1)) +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = plot_ratio_5xfad_long_ba, filename = "Log_ratio_5xfad_gf_long_ba.svg", device = "svg", dpi = "retina", width = 3.5, height = 2)

# linear mixed effect model to test separation
model_ba <- lmer(Ratio ~ host_age * Strain + (1 | host_subject_id), data = data_5xfad_vip_gf_long_ba)
summary(model_ba)


#############
# Export table for Joint-RPCA
#############
data_3xtg_long_join <- data_long_3xtg %>%
  dplyr::left_join(metadata_metabolomics %>% dplyr::select(SampleID, sample_name)) %>%
  dplyr::relocate(sample_name, .after=SampleID) %>% 
  dplyr::mutate(sample_name = gsub("\\.[^.]*$", "", sample_name)) %>% dplyr::select(-SampleID)

#write_csv(x = data_3xtg_long_join, file = "fecal_gf_long_3xtg_metabolomics_join.csv")

data_5xfad_long_join <- data_long_5xfad %>%
  dplyr::left_join(metadata_metabolomics %>% dplyr::select(SampleID, sample_name)) %>%
  dplyr::relocate(sample_name, .after=SampleID) %>% 
  dplyr::mutate(sample_name = gsub("\\.[^.]*$", "", sample_name)) %>% dplyr::select(-SampleID)

#write_csv(x = data_5xfad_long_join, file = "fecal_gf_long_5xfad_metabolomics_join.csv")
