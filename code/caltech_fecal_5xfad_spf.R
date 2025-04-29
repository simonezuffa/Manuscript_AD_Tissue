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
feature_table <- read_csv("data/mice_caltech/fecal_5xfad_spf/gnps_quant_5xfad_spf.csv")
metadata <- read_csv("data/mice_caltech/fecal_5xfad_spf/meta_5xfad_spf.csv")
metadata$tube_id <- as.character(metadata$tube_id)
sample_order <- read.csv("data/mice_caltech/fecal_5xfad_spf/sequence_5xfad_spf.csv") %>% dplyr::select(1,3,4)
annotations <- read.delim("data/mice_caltech/fecal_5xfad_spf/fbmn/nf_output/library/merged_results_with_gnps.tsv") %>%
  dplyr::filter(!str_detect(pattern = "REFRAME", LibraryName)) # remove REFRAME library
annotations$X.Scan. <- as.character(annotations$X.Scan.)
canopus <- read_tsv("data/mice_caltech/fecal_5xfad_spf/canopus_5xfad_spf.tsv") %>% dplyr::select(1,3,5,7,9)
canopus$id <- gsub("^.*sirius_", "", canopus$id)

sample_order$Plate <- gsub("^(.*?):.*", "\\1", sample_order$Position)

info_feature <- feature_table %>% dplyr::select(1:3,7)
colnames(info_feature) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature$Feature <- as.character(info_feature$Feature)

info_feature_complete <- info_feature %>% 
  left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  dplyr::select(1:5,18,24) %>% 
  left_join(canopus %>% distinct(id, .keep_all = TRUE), by = c("Feature" = "id"))

#write_csv(x = info_feature_complete, file = "5xfad_spf_fecal_annotations.csv")

# Extra bile acids libraries and rev cosine
ba_ipsita <- read_csv("data/mice_caltech/fecal_5xfad_spf/merged_Bile_acid_classic_networking.csv") %>% 
  dplyr::filter(SpectrumID %in% annotations$SpectrumID)
ba_rev <- read_tsv("data/mice_caltech/fecal_5xfad_spf/candidate_BA_with_best_annotations.tsv") %>% 
  dplyr::filter(query_id %in% annotations$SpectrumID)


# Data table
data <- feature_table %>% 
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE) %>%
  dplyr::filter(!(str_detect(pattern = "IS_in|IS_fi|IS_en", SampleID)))

data$SampleID <- gsub(".mzML Peak area", "", data$SampleID)

# Metadata
metadata_metabolomics <- data.frame(SampleID = data$SampleID) %>% 
  left_join(metadata, by = c("SampleID" = "tube_id")) %>%
  dplyr::select(1, 3:5, 17, 26, 30, 33, 39, 40) %>% 
  left_join(sample_order)

# Investigate total peak area
data_TIC <- data.frame(TIC = rowSums(data %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

data_TIC %>% dplyr::filter(mouse_strain %in% c("WT 5XFAD", "Het 5XFAD")) %>%
  ggscatter("Run_Order", "TIC", add = "reg.line") + ylim(0, 9e8) +
  stat_cor()

# Check sample type
sample_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "36", SampleID)) %>% summarise(mean(TIC))
pool_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Pool", SampleID)) %>% summarise(mean(TIC))
six_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "QCmix", SampleID)) %>% summarise(mean(TIC))
blank_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "BLANK", SampleID)) %>% summarise(mean(TIC))
SRM_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "SRM", SampleID)) %>% summarise(mean(TIC))

ratio_tic_pb <- pool_tic/sample_tic
ratio_tic_sb <- sample_tic/blank_tic
ratio_tic_6p <- pool_tic/six_tic

# TIC of QCpool is ~3 times the one of samples --> probably error in sample preparation 
# Divide intensity by 2.5 but not for the internal standard cause it was reconstituted in one (checked)

data_nopool <- data %>% dplyr::filter(!(str_detect(pattern = "Pool", SampleID)))
pool_data <- data %>% dplyr::filter(str_detect(pattern = "Pool", SampleID))

pool_data_fix <- pool_data %>% column_to_rownames("SampleID") %>% 
  dplyr::mutate(across(everything(), ~./2.5)) %>%
  dplyr::mutate(`5305` = pool_data$`5305`) %>% rownames_to_column("SampleID")

data_fix <- rbind(data_nopool, pool_data_fix)

# Calculate TIC
data_TIC <- data.frame(TIC = rowSums(data_fix %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

data_TIC %>% ggscatter("Run_Order", "TIC", add = "reg.line") +
  stat_cor() + ylim(0, 9e8)

# Now it is fixed. One problem is that the intensity of the samples is not very high compared to blanks

# Check sample type
sample_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "36", SampleID)) %>% summarise(mean(TIC))
pool_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Pool", SampleID)) %>% summarise(mean(TIC))
six_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "QCmix", SampleID)) %>% summarise(mean(TIC))
blank_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "BLANK", SampleID)) %>% summarise(mean(TIC))
SRM_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "SRM", SampleID)) %>% summarise(mean(TIC))

ratio_tic_pb <- pool_tic/sample_tic
ratio_tic_sb <- sample_tic/blank_tic
ratio_tic_6p <- pool_tic/six_tic

# Check TIC overtime for QCpool, QCmix, Blank, SRM
data_TIC %>% dplyr::filter(str_detect(pattern = "Pool", SampleID)) %>% 
  ggscatter("Run_Order", "TIC", add = "reg.line") + ylim(0, 4e8) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "QCmix", SampleID)) %>% 
  ggscatter("Run_Order", "TIC", add = "reg.line") + ylim(0, 6e8) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "BLANK", SampleID)) %>% 
  ggscatter("Run_Order", "TIC", add = "reg.line") + ylim(0, 4e8) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "SRM", SampleID)) %>% 
  ggscatter("Run_Order", "TIC", add = "reg.line") + ylim(0, 4e8) +
  stat_cor()
# Blanks and QCpool TICs are declining a bit through the run --> ok

# Check internal standard and remove sample where there was a shift
fbmn_IS <- annotations %>% dplyr::filter(str_detect(Compound_Name, regex("sulf", ignore_case = TRUE))) %>% 
  distinct(X.Scan., .keep_all = TRUE) %>% dplyr::filter(Organism != "BILELIB19")

# Extract IS features for the table
table_IS <- data_fix %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% 
  dplyr::filter(ID %in% fbmn_IS$X.Scan.) %>% column_to_rownames("ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SampleID") %>% dplyr::filter(!(str_detect(SampleID, "SRM|QCmix|Pool|BLANK|SRM|IS"))) %>%
  dplyr::select(SampleID, `5305`) %>% left_join(metadata_metabolomics) %>% 
  dplyr::filter(`5305` < 3000000) %>% dplyr::filter(`5305` > 1200000) # remove samples with poor IS acquisition

sample_is_to_remove <- data_fix %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% 
  dplyr::filter(ID %in% fbmn_IS$X.Scan.) %>% column_to_rownames("ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SampleID") %>% dplyr::filter(!(str_detect(SampleID, "SRM|QCmix|Pool|BLANK|SRM|IS"))) %>%
  dplyr::select(SampleID, `5305`) %>% left_join(metadata_metabolomics) %>% 
  dplyr::filter(`5305` > 3000000 | `5305` < 1200000)

colnames(table_IS)[2] <- "Sulfadimethoxine"

table_IS %>% 
  ggscatter(x = "Run_Order", y = "Sulfadimethoxine", add = c("reg.line")) + ylim(0, 4e6) +
  stat_cor() 

table_IS %>% 
  ggbarplot(x = "Run_Order", y = "Sulfadimethoxine", xlab = "Run Order", 
            ylab = "Peak Area Sulfadimethoxine", title = "Internal Standard Acquisition") +
  geom_hline(yintercept = mean(table_IS$Sulfadimethoxine, na.rm = TRUE), linetype = "dashed", color = "blue")

cv_is <- sd(table_IS$Sulfadimethoxine)/mean(table_IS$Sulfadimethoxine)


# Check features in blanks and pools
data_blank <- data_fix %>% dplyr::filter(str_detect(pattern = "BLANK", SampleID))
data_pools <- data_fix %>% dplyr::filter(str_detect(pattern = "Pool", SampleID))
data_sixmix <- data_fix %>% dplyr::filter(str_detect(pattern = "QCmix", SampleID))
data_srm <- data_fix %>% dplyr::filter(str_detect(pattern = "SRM", SampleID))

# Blanks
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

# Pool
pools_feature_info <- data.frame(Feature = colnames(data_pools)[-1],
                                 Mean_pool = data_pools %>% column_to_rownames("SampleID") %>% colMeans(), 
                                 SD_pool =  data_pools %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
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
data_clean <- data_fix %>% dplyr::select(-c(feature_to_remove$Feature)) %>%
  dplyr::filter(!(SampleID %in% c("363243385", "363243387", "363141197"))) %>%
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")

# Features to be removed Pool/QCmix < 5
feature_to_remove_mix <- sixmix_feature_info %>% left_join(pools_feature_info) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% 
  dplyr::mutate(Pool_Mix = Mean_pool/Mean_sixmix) %>% 
  dplyr::filter(Pool_Mix < 5 | is.na(Pool_Mix)) %>% 
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))

# Data with blank2 removal
data_clean2 <- data_clean %>% dplyr::select(-c(feature_to_remove_mix$Feature))

# Remove feature before 0.2 minutes and after 8 minutes
feature_to_remove_rt <- info_feature_complete %>% dplyr::filter(RT < 0.2 | RT > 8) %>%
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))  %>%
  dplyr::filter(!(Feature %in% feature_to_remove_mix$Feature)) %>%
  dplyr::filter(Feature != "14823")

# Final cleaned table
data_clean3 <- data_clean2 %>% dplyr::select(-c(feature_to_remove_rt$Feature)) %>%
  dplyr::filter(!(SampleID %in% sample_is_to_remove$SampleID))


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
  dplyr::filter(!(str_detect(pattern = "mix|BLANK|SRM|IS_|Pool", SampleID)))

# PCA raw data
PCA_raw <- mixOmics::pca(data_sample %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
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

# RCLR transformation
data_sample_clr <- decostand(data_sample %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

i <- "intervention"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = "Fecal - SPF 5xFAD", palette = c("#B7E6A5", "#007188"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = PCA_plot, filename = "5xfad_spf_fecal_pca_metabolome.svg", device = "svg", dpi = "retina", height = 2, width = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_sample_clr, method = "euclidean")
disper_intervention <- betadisper(dist_metabolites, PCA_whole_scores$intervention)
anova(disper_intervention)
permanova <- adonis2(dist_metabolites ~ intervention, PCA_whole_scores, na.action = na.omit)

# Work separately on the two different studies cause there is a strong clustering 


# Check sacrifice and longitudinal data separately
sample_sac <- metadata_metabolomics %>% dplyr::filter(intervention == "Sacrifice")
sample_long <- metadata_metabolomics %>% dplyr::filter(intervention != "Sacrifice")

sample_wt <- metadata_metabolomics %>% dplyr::filter(mouse_strain == "WT 5XFAD")
sample_mut <- metadata_metabolomics %>% dplyr::filter(mouse_strain == "Het 5XFAD")

sample_male <- metadata_metabolomics %>% dplyr::filter(sex == "male")
sample_female <- metadata_metabolomics %>% dplyr::filter(sex == "female")


##################
# SACRIFICE DATA #
##################
data_sacrifice <- data_sample %>% 
  dplyr::filter(SampleID %in% sample_sac$SampleID)

# RCLR transformation
data_sacrifice_clr <- decostand(data_sacrifice %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sacrifice_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_sac_plots <- list()

for (i in c("cohort", "host_age", "mouse_strain", "sex", "strain_sex", "Plate")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = "PCA - Sacrifice",
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_sac_plots[[i]] <- PCA_plot
  
}

PCA_sac_plots_final <- wrap_plots(PCA_sac_plots, nrow = 3)

i <- "mouse_strain"

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

#ggsave(plot = PCA_plot, filename = "5xfad_spf_fecal_sacrifice_pca_metabolome.svg", device = "svg", dpi = "retina", height = 2, width = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_sacrifice_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$mouse_strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ mouse_strain * sex * host_age + 
                       cohort + Plate, PCA_whole_scores, na.action = na.omit)

#write_csv(x = data_sacrifice, file = "fecal_spf_sacrifice_5xfad_metabolome.csv")
#write_csv(x = data_sacrifice_clr, file = "fecal_spf_sacrifice_5xfad_metabolome_rclr.csv")
#write_csv(x = metadata_metabolomics %>% dplyr::filter(SampleID %in% data_sacrifice$SampleID), file = "fecal_spf_sacrifice_5xfad_metadata.csv")


# PLSDA 5xFAD - Strain
PLSDA_strain <- mixOmics::plsda(data_sacrifice_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                PCA_whole_scores$mouse_strain, ncomp = 3, scale = TRUE)
PLSDA_strain_scores <- data.frame(PLSDA_strain$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_strain_plot <- PLSDA_strain_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "mouse_strain", alpha = 0.6, title = "PLSDA - Fecal 5xFAD SPF - Strain",
            xlab = paste("Component 1 (", round(PLSDA_strain$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_strain$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_strain_scores %>% group_by(mouse_strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_strain <- plotLoadings(PLSDA_strain, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_strain <- perf(PLSDA_strain, validation = "Mfold", folds = 4, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_strain, legend = FALSE)

VIPs_strain <- as.data.frame(mixOmics::vip(PLSDA_strain))
VIPs_strain_filter <- dplyr::filter(VIPs_strain, VIPs_strain$comp1 > 1)
VIPs_strain_filter$ID <- rownames(VIPs_strain_filter)
VIPs_strain_select <- VIPs_strain_filter %>% dplyr::select(ID, comp1)
VIPs_strain_Load <- VIPs_strain_select %>% 
  left_join(Loadings_strain, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

# PLSDA 5xFAD - Sex
PLSDA_sex <- mixOmics::plsda(data_sacrifice_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                             PCA_whole_scores$sex, ncomp = 3, scale = TRUE)
PLSDA_sex_scores <- data.frame(PLSDA_sex$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_sex_plot <- PLSDA_sex_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "sex", alpha = 0.6, title = "PLSDA - Fecal 5xFAD SPF - Sex",
            xlab = paste("Component 1 (", round(PLSDA_sex$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_sex$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_sex_scores %>% group_by(sex) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = sex), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_sex <- plotLoadings(PLSDA_sex, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_sex <- perf(PLSDA_sex, validation = "Mfold", folds = 4, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_sex, legend = FALSE)

VIPs_sex <- as.data.frame(mixOmics::vip(PLSDA_sex))
VIPs_sex_filter <- dplyr::filter(VIPs_sex, VIPs_sex$comp1 > 1)
VIPs_sex_filter$ID <- rownames(VIPs_sex_filter)
VIPs_sex_select <- VIPs_sex_filter %>% dplyr::select(ID, comp1)
VIPs_sex_Load <- VIPs_sex_select %>% 
  left_join(Loadings_sex, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))


# Remove feature influenced by sex
VIPs_strain_5xfad_no_sex <- VIPs_strain_Load %>% 
  anti_join(VIPs_sex_Load %>% dplyr::select(ID))

#write_csv(x = VIPs_strain_5xfad_no_sex, "5xfad_spf_fecal_plsda_vip_all_samples_genotype.csv")


# Stratify by genotype and sex and look at differences if
data_sac_het <- data_sacrifice %>% 
  dplyr::filter(SampleID %in% sample_mut$SampleID)
data_sac_wt <- data_sacrifice %>% 
  dplyr::filter(SampleID %in% sample_wt$SampleID)
data_sac_male <- data_sacrifice %>% 
  dplyr::filter(SampleID %in% sample_male$SampleID)
data_sac_female <- data_sacrifice %>% 
  dplyr::filter(SampleID %in% sample_female$SampleID)

# RCLR transformation
data_sac_het_clr <- decostand(data_sac_het %>% column_to_rownames("SampleID"), method = "rclr")
data_sac_wt_clr <- decostand(data_sac_wt %>% column_to_rownames("SampleID"), method = "rclr")
data_sac_male_clr <- decostand(data_sac_male %>% column_to_rownames("SampleID"), method = "rclr")
data_sac_female_clr <- decostand(data_sac_female %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_sac_het <- mixOmics::pca(data_sac_het_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                             ncomp = 2, center = TRUE, scale = TRUE)
PCA_sac_wt <- mixOmics::pca(data_sac_wt_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                            ncomp = 2, center = TRUE, scale = TRUE)
PCA_sac_male <- mixOmics::pca(data_sac_male_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                              ncomp = 2, center = TRUE, scale = TRUE)
PCA_sac_female <- mixOmics::pca(data_sac_female_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                ncomp = 2, center = TRUE, scale = TRUE)

PCA_sac_het_scores <- data.frame(PCA_sac_het$variates$X) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)
PCA_sac_wt_scores <- data.frame(PCA_sac_wt$variates$X) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)
PCA_sac_male_scores <- data.frame(PCA_sac_male$variates$X) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)
PCA_sac_female_scores <- data.frame(PCA_sac_female$variates$X) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_sac_het_plot <- PCA_sac_het_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "sex", alpha = 0.6,
            title = "PCA - Het Sacrifice",
            xlab = paste("PC1 (", round(PCA_sac_het$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sac_het$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sac_het_scores %>% group_by(sex) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2)), 
             aes(x = PC1, y = PC2, color = sex), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), 
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + 
  coord_fixed()

PCA_sac_wt_plot <- PCA_sac_wt_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "sex", alpha = 0.6,
            title = "PCA - WT Sacrifice",
            xlab = paste("PC1 (", round(PCA_sac_wt$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sac_wt$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sac_wt_scores %>% group_by(sex) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2)), 
             aes(x = PC1, y = PC2, color = sex), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), 
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + 
  coord_fixed()

PCA_sac_male_plot <- PCA_sac_male_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "mouse_strain", alpha = 0.6,
            title = "PCA - Male Sacrifice",
            xlab = paste("PC1 (", round(PCA_sac_male$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sac_male$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sac_male_scores %>% group_by(mouse_strain) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2)), 
             aes(x = PC1, y = PC2, color = mouse_strain), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), 
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + 
  coord_fixed()

PCA_sac_female_plot <- PCA_sac_female_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "mouse_strain", alpha = 0.6,
            title = "PCA - Female Sacrifice",
            xlab = paste("PC1 (", round(PCA_sac_female$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sac_female$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sac_female_scores %>% group_by(mouse_strain) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2)), 
             aes(x = PC1, y = PC2, color = mouse_strain), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), 
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + 
  coord_fixed()

PCA_sac_subplots <- wrap_plots(PCA_sac_het_plot, PCA_sac_wt_plot, PCA_sac_male_plot, PCA_sac_female_plot, nrow = 2)

# PLSDA
PLSDA_sac_het <- mixOmics::plsda(data_sac_het_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                 PCA_sac_het_scores$sex, ncomp = 2, scale = TRUE)
PLSDA_sac_wt <- mixOmics::plsda(data_sac_wt_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                PCA_sac_wt_scores$sex, ncomp = 3, scale = TRUE)
PLSDA_sac_male <- mixOmics::plsda(data_sac_male_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                  PCA_sac_male_scores$mouse_strain, ncomp = 3, scale = TRUE)
PLSDA_sac_female <- mixOmics::plsda(data_sac_female_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                    PCA_sac_female_scores$mouse_strain, ncomp = 2, scale = TRUE)

PLSDA_sac_het_scores <- data.frame(PLSDA_sac_het$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)
PLSDA_sac_wt_scores <- data.frame(PLSDA_sac_wt$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)
PLSDA_sac_male_scores <- data.frame(PLSDA_sac_male$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)
PLSDA_sac_female_scores <- data.frame(PLSDA_sac_female$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_sac_het_plot <- PLSDA_sac_het_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "sex", alpha = 0.6, title = "PLSDA - Het Sacrifice Sex",
            xlab = paste("Component 1 (", round(PLSDA_sac_het$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_sac_het$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_sac_het_scores %>% group_by(sex) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = sex), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_sac_het <- plotLoadings(PLSDA_sac_het, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_sac_het <- perf(PLSDA_sac_het, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_sac_het, legend = FALSE)

PLSDA_sac_wt_plot <- PLSDA_sac_wt_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "sex", alpha = 0.6, title = "PLSDA - WT Sacrifice Sex",
            xlab = paste("Component 1 (", round(PLSDA_sac_wt$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_sac_wt$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_sac_wt_scores %>% group_by(sex) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = sex), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_sac_wt <- plotLoadings(PLSDA_sac_wt, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_sac_wt <- perf(PLSDA_sac_wt, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_sac_wt, legend = FALSE)

PLSDA_sac_male_plot <- PLSDA_sac_male_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "mouse_strain", alpha = 0.6, title = "PLSDA - Sacrifice Male Genotype",
            xlab = paste("Component 1 (", round(PLSDA_sac_male$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_sac_male$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_sac_male_scores %>% group_by(mouse_strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_sac_male <- plotLoadings(PLSDA_sac_male, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_sac_male <- perf(PLSDA_sac_male, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_sac_male, legend = FALSE)

PLSDA_sac_female_plot <- PLSDA_sac_female_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "mouse_strain", alpha = 0.6, title = "PLSDA - Sacrifice Female Genotype",
            xlab = paste("Component 1 (", round(PLSDA_sac_female$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_sac_female$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_sac_female_scores %>% group_by(mouse_strain) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_strain), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_sac_female <- plotLoadings(PLSDA_sac_female, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_sac_female <- perf(PLSDA_sac_female, validation = "Mfold", folds = 5, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_sac_female, legend = FALSE)


VIPs_sac_het <- as.data.frame(mixOmics::vip(PLSDA_sac_het))
VIPs_sac_het_filter <- dplyr::filter(VIPs_sac_het, VIPs_sac_het$comp1 > 1)
VIPs_sac_het_filter$ID <- rownames(VIPs_sac_het_filter)
VIPs_sac_het_select <- VIPs_sac_het_filter %>% dplyr::select(ID, comp1)
VIPs_sac_het_Load <- VIPs_sac_het_select %>% 
  left_join(Loadings_sac_het, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

VIPs_sac_wt <- as.data.frame(mixOmics::vip(PLSDA_sac_wt))
VIPs_sac_wt_filter <- dplyr::filter(VIPs_sac_wt, VIPs_sac_wt$comp1 > 1)
VIPs_sac_wt_filter$ID <- rownames(VIPs_sac_wt_filter)
VIPs_sac_wt_select <- VIPs_sac_wt_filter %>% dplyr::select(ID, comp1)
VIPs_sac_wt_Load <- VIPs_sac_wt_select %>% 
  left_join(Loadings_sac_wt, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

VIPs_sac_male <- as.data.frame(mixOmics::vip(PLSDA_sac_male))
VIPs_sac_male_filter <- dplyr::filter(VIPs_sac_male, VIPs_sac_male$comp1 > 1)
VIPs_sac_male_filter$ID <- rownames(VIPs_sac_male_filter)
VIPs_sac_male_select <- VIPs_sac_male_filter %>% dplyr::select(ID, comp1)
VIPs_sac_male_Load <- VIPs_sac_male_select %>% 
  left_join(Loadings_sac_male, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

VIPs_sac_female <- as.data.frame(mixOmics::vip(PLSDA_sac_female))
VIPs_sac_female_filter <- dplyr::filter(VIPs_sac_female, VIPs_sac_female$comp1 > 1)
VIPs_sac_female_filter$ID <- rownames(VIPs_sac_female_filter)
VIPs_sac_female_select <- VIPs_sac_female_filter %>% dplyr::select(ID, comp1)
VIPs_sac_female_Load <- VIPs_sac_female_select %>% 
  left_join(Loadings_sac_female, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))


# Sample matching 
list_5xfad_genotype <- list(
  `Male 5xFAD` = (VIPs_sac_male_Load %>% dplyr::filter(GroupContrib == "Het 5XFAD"))$ID,
  `Male WT` =  (VIPs_sac_male_Load %>% dplyr::filter(GroupContrib == "WT 5XFAD"))$ID,
  `Female 5xFAD` = (VIPs_sac_female_Load %>% dplyr::filter(GroupContrib == "Het 5XFAD"))$ID,
  `Female WT` = (VIPs_sac_female_Load %>% dplyr::filter(GroupContrib == "WT 5XFAD"))$ID)

genotype_5xfad <- UpSetR::upset(fromList(list_5xfad_genotype), nsets = 4, nintersects = 7, 
                                point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                                sets = c("Male 5xFAD", "Male WT", "Female 5xFAD", "Female WT"),
                                queries = list(list(query = intersects, params = list("Female 5xFAD", 
                                                                                      "Male 5xFAD"), color = "steelblue1", active = T), 
                                               list(query = intersects, params = list("Female WT", 
                                                                                      "Male WT"),color = "orange", active = T)))

#pdf(file="5xfad_spf_fecal_upset_genotype.pdf", width = 3, height = 2.5)
#genotype_5xfad
#dev.off()

list_5xfad_sex <- list(
  `5xFAD Male` = (VIPs_sac_het_Load %>% dplyr::filter(GroupContrib == "male"))$ID,
  `5xFAD Female` =  (VIPs_sac_het_Load %>% dplyr::filter(GroupContrib == "female"))$ID,
  `WT Male` = (VIPs_sac_wt_Load %>% dplyr::filter(GroupContrib == "male"))$ID,
  `WT Female` = (VIPs_sac_wt_Load %>% dplyr::filter(GroupContrib == "female"))$ID)

sex_5xfad <- UpSetR::upset(fromList(list_5xfad_sex), nsets = 4, nintersects = 6, 
                                point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                                sets = c("5xFAD Male", "5xFAD Female", "WT Male", "WT Female"),
                                queries = list(list(query = intersects, params = list("WT Male", 
                                                                                      "5xFAD Male"), color = "steelblue1", active = T), 
                                               list(query = intersects, params = list("5xFAD Female", 
                                                                                      "WT Female"),color = "orange", active = T)))

#pdf(file="5xfad_spf_fecal_upset_sex.pdf", width = 3, height = 2.5)
#sex_5xfad
#dev.off()

list_5xfad_combined <- list(
  `Male 5xFAD` = (VIPs_sac_male_Load %>% dplyr::filter(GroupContrib == "Het 5XFAD"))$ID,
  `Male WT` =  (VIPs_sac_male_Load %>% dplyr::filter(GroupContrib == "WT 5XFAD"))$ID,
  `Female 5xFAD` = (VIPs_sac_female_Load %>% dplyr::filter(GroupContrib == "Het 5XFAD"))$ID,
  `Female WT` = (VIPs_sac_female_Load %>% dplyr::filter(GroupContrib == "WT 5XFAD"))$ID,
  `5xFAD Male` = (VIPs_sac_het_Load %>% dplyr::filter(GroupContrib == "male"))$ID,
  `5xFAD Female` =  (VIPs_sac_het_Load %>% dplyr::filter(GroupContrib == "female"))$ID,
  `WT Male` = (VIPs_sac_wt_Load %>% dplyr::filter(GroupContrib == "male"))$ID,
  `WT Female` = (VIPs_sac_wt_Load %>% dplyr::filter(GroupContrib == "female"))$ID)

combined_5xfad <- UpSetR::upset(fromList(list_5xfad_combined), nsets = 8, nintersects = NA,
                               point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                               sets = c("5xFAD Male", "5xFAD Female", "WT Male", "WT Female", "Male 5xFAD", "Male WT", "Female 5xFAD", "Female WT"),
                               queries = list(list(query = intersects, params = list("WT Male", 
                                                                                     "5xFAD Male"), color = "orange", active = T), 
                                              list(query = intersects, params = list("5xFAD Female", 
                                                                                     "WT Female"),color = "orange", active = T), 
                                              list(query = intersects, params = list("Female 5xFAD", 
                                                                                     "Male 5xFAD"), color = "steelblue1", active = T), 
                                              list(query = intersects, params = list("Female WT", 
                                                                                     "Male WT"),color = "steelblue1", active = T)))

#pdf(file="5xfad_spf_fecal_upset_genotype_sex.pdf", width = 3, height = 2.5)
#combined_5xfad
#dev.off()


# Check features of interest
vip_5xfad_spf_wt <- VIPs_strain_Load %>% dplyr::filter(GroupContrib == "WT 5XFAD")
vip_5xfad_spf_mut <- VIPs_strain_Load %>% dplyr::filter(GroupContrib == "Het 5XFAD")

data_5xfad_vip_spf <- data_sacrifice %>%
  dplyr::select("SampleID", VIPs_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_5xfad_spf_wt$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_5xfad_spf_mut$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_metabolomics)

plot_ratio_feces_5xfad <- data_5xfad_vip_spf %>% 
  dplyr::mutate(mouse_strain = factor(mouse_strain, levels = c("WT 5XFAD", "Het 5XFAD"))) %>%
  ggboxplot(y = "Ratio", x = "mouse_strain", add = "jitter", ylab = "Ln(WT/Mut)",
            add.params = list(color = "mouse_strain", alpha = 0.6), legend = "none",
            palette = c("#E69A8DFF", "#5F4B8BFF"), xlab = "Feces",
            title = "Differential features from PLS-DA models") + 
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = plot_ratio_feces_5xfad, filename = "5xfad_feces_spf_ratio.svg", device = "svg", dpi = "retina", width = 1.5, height = 2.5)
#write_csv(x = data_5xfad_vip_spf, file = "feces_spf_ratio_5xfad.csv")

# Cliff's Delta
cliff.delta((data_5xfad_vip_spf %>%
               dplyr::filter(mouse_strain == "WT 5XFAD"))$Ratio, 
            (data_5xfad_vip_spf %>%
               dplyr::filter(mouse_strain == "Het 5XFAD"))$Ratio)


# Extract features of interest for MN
fecal_spf_5xfad_interest <- VIPs_strain_Load %>% 
  dplyr::select(ID, comp1, GroupContrib) %>% 
  left_join(VIPs_sac_male_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>% 
  left_join(VIPs_sac_female_Load %>% 
              dplyr::select(ID, comp1, GroupContrib), by = c("ID" = "ID")) %>%
  dplyr::mutate(VIP_mean = rowMeans(dplyr::select(., comp1.x, comp1.y, comp1), na.rm = TRUE))

colnames(fecal_spf_5xfad_interest) <- c("ID", "VIP1", "All_genotype", "VIP2", "Male_genotype",
                                       "VIP3", "Female_genotype", "VIP_mean")

fecal_spf_5xfad_interest <- fecal_spf_5xfad_interest %>%   
  dplyr::mutate(Genotype = case_when(Male_genotype == Female_genotype ~ "Yes",
                                     TRUE ~ "No")) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature"))

#write_csv(x = fecal_spf_5xfad_interest, file = "5xfad_spf_fecal_features_info_mn.csv")


##########
# Export table for Joint-RPCA
##########
data_fecal_join <- data_sacrifice %>%
  dplyr::left_join(metadata_metabolomics %>% dplyr::select(SampleID, anonymized_name)) %>%
  dplyr::relocate(anonymized_name, .after=SampleID) %>% 
  dplyr::mutate(anonymized_name = gsub("\\.[^.]*$", "", anonymized_name)) %>% dplyr::select(-SampleID)

#write_csv(x = data_fecal_join, file = "fecal_spf_sac_5xfad_metabolomics_join.csv")


#####################
# LONGITUDINAL DATA #
#####################
data_longitudinal <- data_sample %>% 
  dplyr::filter(SampleID %in% sample_long$SampleID)

metadata_longitudinal <- metadata_metabolomics %>% 
  dplyr::filter(intervention == "Longitudinal") %>%
  dplyr::mutate_at("host_age", as.numeric) %>% arrange(host_age)

# Plot longitudinal metadata information
meta_long_5xfad_plot <- metadata_longitudinal %>% 
  group_by(host_age, mouse_strain) %>%   
  summarise(count = n(), .groups = "drop") %>%
  ggplot(aes(x = factor(host_age), y = count, fill = mouse_strain)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#5F4B8BFF", "#E69A8DFF")) +
  labs(x = "Age (months)", y = "Count", fill = "Genotype",
       title = "Fecal samples in 5xFAD SPF") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = meta_long_5xfad_plot, filename = "Metadata_fecal_5xfad_long.svg", device = "svg", dpi = "retina", width = 2.5, height = 2)


# RCLR transformation
data_longitudinal_clr <- decostand(data_longitudinal %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_longitudinal_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)
PCA_whole_scores$host_age <- as.numeric(PCA_whole_scores$host_age)

PCA_long_plots <- list()

for (i in c("cohort", "host_age", "host_subject_id", "mouse_strain", "sex", "strain_sex", "Plate")) {
  
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
dist_metabolites <- vegdist(data_longitudinal_clr, method = "euclidean")
disper_genotype <- betadisper(dist_metabolites, PCA_whole_scores$mouse_strain)
anova(disper_genotype)
permanova <- adonis2(dist_metabolites ~ mouse_strain + sex + host_age + host_subject_id, PCA_whole_scores, na.action = na.omit)

i <- "mouse_strain"

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


# Extract features and generate log ratio which was established in the Sacrifice cohort
data_5xfad_vip_spf_long <- data_longitudinal %>%
  dplyr::select("SampleID", VIPs_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_5xfad_spf_wt$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_5xfad_spf_mut$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_longitudinal)

plot_ratio_5xfad_long <- data_5xfad_vip_spf_long %>% 
  dplyr::mutate_at("host_age", as.numeric) %>%
  ggscatter(x = "host_age", y = "Ratio", color = "mouse_strain", add = "loess", ylab = "Ln(TopWT/TopMut)",
            xlab = "Age (months)", legend = "none", palette = c("#5F4B8BFF", "#E69A8DFF"), alpha = 0.4,
            title = "Natural Log Ratio Across Time in 5xFAD SPF") + 
  scale_x_continuous(breaks = seq(min(data_5xfad_vip_spf_long$host_age, na.rm = TRUE),
                                  max(data_5xfad_vip_spf_long$host_age, na.rm = TRUE), by = 1)) +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = plot_ratio_5xfad_long, filename = "Log_ratio_5xfad_spf_long.svg", device = "svg", dpi = "retina", width = 4, height = 2)

# linear mixed effect model to test separation
library(lmerTest)
model <- lmer(Ratio ~ host_age * mouse_strain + (1 | host_subject_id), data = data_5xfad_vip_spf_long)
summary(model)


# check bile acid ratio only
ba_massql <- read_tsv("data/mice_caltech/fecal_5xfad_spf/ba_massql_5xfad_feces.tsv")
ba_massql_no <- ba_massql%>%
  dplyr::filter(str_detect(pattern = "Did not pass", query_validation))

vip_5xfad_spf_mut_ba <- vip_5xfad_spf_mut %>% 
  dplyr::filter(str_detect(pattern = "butanoic acid|pentanoic acid|bile acid|dihydroxyoxane-2-carboxylic acid", Compound_Name)) %>%
  dplyr::filter(!(ID %in% ba_massql_no$`#Scan#`))
vip_5xfad_spf_wt_ba <- vip_5xfad_spf_wt %>% 
  dplyr::filter(str_detect(pattern = "butanoic acid|pentanoic acid|bile acid|dihydroxyoxane-2-carboxylic acid", Compound_Name)) %>%
  dplyr::filter(!(ID %in% ba_massql_no$`#Scan#`))

data_5xfad_vip_spf_long_ba <- data_longitudinal %>%
  dplyr::select("SampleID", VIPs_strain_Load$ID) %>%
  dplyr::mutate(WT = rowSums(select(., vip_5xfad_spf_wt_ba$ID))) %>%
  dplyr::mutate(Mut = rowSums(select(., vip_5xfad_spf_mut_ba$ID))) %>%
  dplyr::mutate(Ratio = log(WT/Mut)) %>%
  left_join(metadata_longitudinal)

plot_ratio_5xfad_long_ba <- data_5xfad_vip_spf_long_ba %>% 
  dplyr::mutate_at("host_age", as.numeric) %>%
  ggscatter(x = "host_age", y = "Ratio", color = "mouse_strain", add = "loess", ylab = "Ln(WTTop/MutTop)",
            xlab = "Age (months)", legend = "none", palette = c("#5F4B8BFF", "#E69A8DFF"), alpha = 0.4,
            title = "Bile Acid Ratio in Fecal 5xfad SPF") + 
  scale_x_continuous(breaks = seq(min(data_5xfad_vip_spf_long_ba$host_age, na.rm = TRUE),
                                  max(data_5xfad_vip_spf_long_ba$host_age, na.rm = TRUE), by = 1)) +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = plot_ratio_5xfad_long_ba, filename = "Log_ratio_5xfad_spf_long_ba.svg", device = "svg", dpi = "retina", width = 3.5, height = 2)

# linear mixed effect model to test separation
model_ba <- lmer(Ratio ~ host_age * mouse_strain + (1 | host_subject_id), data = data_5xfad_vip_spf_long_ba)
summary(model_ba)


##########
# Export table for Joint-RPCA
##########
data_fecal_long_join <- data_longitudinal %>%
  dplyr::left_join(metadata_metabolomics %>% dplyr::select(SampleID, anonymized_name)) %>%
  dplyr::relocate(anonymized_name, .after=SampleID) %>% 
  dplyr::mutate(anonymized_name = gsub("\\.[^.]*$", "", anonymized_name)) %>% dplyr::select(-SampleID)

#write_csv(x = data_fecal_long_join, file = "fecal_spf_long_5xfad_metabolomics_join.csv")
