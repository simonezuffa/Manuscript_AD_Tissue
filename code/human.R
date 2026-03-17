setwd("~/Desktop/Manuscript_AD_Tissue")

library(tidyverse)
library(mixOmics)
library(ggpubr)
library(vegan)
library(caret)
library(limma)
library(patchwork)
library(rstatix)
library(effsize)
library(readxl)

# Read data
adrc <- read_csv("data/human/adrc_plasma/mzmine/adrc_quant.csv")
t1000 <- read_csv("data/human/t1000_serum/mzmine/t1000_quant.csv")
mars <- read_csv("data/human/mars_serum/mzmine/wisconsin_quant.csv")

info_feature_adrc <- adrc %>% dplyr::select(1:3,7)
colnames(info_feature_adrc) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature_adrc$Feature <- as.character(info_feature_adrc$Feature)

info_feature_t1000 <- t1000 %>% dplyr::select(1:3,7)
colnames(info_feature_t1000) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature_t1000$Feature <- as.character(info_feature_t1000$Feature)

info_feature_mars <- mars %>% dplyr::select(1:3,7)
colnames(info_feature_mars) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature_mars$Feature <- as.character(info_feature_mars$Feature)

# My feature of interest has RT ~ 2.20 min in serum
# adrc        7132
# t1000       3031
# win         2458


# Read metadata - unfortunately I cannot share but can be requested via https://adknowledgeportal.synapse.org/
adrc_meta <- read_csv("data/human/adrc_plasma/metadata/")
adrc_meta$UCSD <- as.character(adrc_meta$UCSD)

t1000_meta <- read_csv("data/human/t1000_serum/metadata/")
t1000_key <- read_excel("data/human/t1000_serum/metadata/x")
t1000_meta_join <- t1000_key %>% left_join(t1000_meta, by = c("Patient ID" = "id"))

mars_meta <- read_csv("data/human/mars_serum/metadata/")
mars_csf_biomarkers <- read.delim("data/human/mars_serum/metadata/")
mars_csf_extra <- read_csv("data/human/mars_serum/metadata/")
mars_key <- read_excel("data/human/mars_serum/metadata/") %>%
  dplyr::select(fecalMarsUniqueID, serumSpecimenBarCode, serumMARSID)
mars_meta_join <- mars_key %>% left_join(mars_meta, by = c("fecalMarsUniqueID" = "uniqueid")) %>%
  left_join(mars_csf_biomarkers, by = c("fecalMarsUniqueID" = "unique_id"))
mars_meta_join$serumSpecimenBarCode <- as.character(mars_meta_join$serumSpecimenBarCode)


# Data tables
data_adrc <- adrc %>%
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE)
data_adrc$SampleID <- gsub(".mzML Peak area", "", data_adrc$SampleID)

data_t1000 <- t1000 %>%
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE)
data_t1000$SampleID <- gsub(".mzXML Peak area", "", data_t1000$SampleID)

data_mars <- mars %>%
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE)
data_mars$SampleID <- gsub(".mzML Peak area", "", data_mars$SampleID)


# Check TIC 
data_adrc_TIC <- data.frame(TIC = rowSums(data_adrc %>% column_to_rownames("SampleID"))) %>%
  rownames_to_column("SampleID") %>% arrange(TIC) %>% 
  dplyr::mutate(Type = case_when(str_detect(pattern = "pool", SampleID) ~ "QC_pool",
                                 str_detect(pattern = "srm", SampleID) ~ "SRM",
                                 str_detect(pattern = "ACN", SampleID) ~ "Blank",
                                 str_detect(pattern = "mix", SampleID) ~ "QC_mix",
                                 str_detect(pattern = "Sol", SampleID) ~ "Solvent",
                                 TRUE ~ "Sample"))

data_t1000_TIC <- data.frame(TIC = rowSums(data_t1000 %>% column_to_rownames("SampleID"))) %>%
  rownames_to_column("SampleID") %>% arrange(TIC) %>% 
  dplyr::mutate(Type = case_when(str_detect(pattern = "Blank", SampleID) ~ "Blank",
                                 str_detect(pattern = "QC", SampleID) ~ "QC",
                                 TRUE ~ "Sample"))

data_mars_TIC <- data.frame(TIC = rowSums(data_mars %>% column_to_rownames("SampleID"))) %>%
  rownames_to_column("SampleID") %>% arrange(TIC) %>% 
  dplyr::mutate(Type = case_when(str_detect(pattern = "Blank", SampleID) ~ "Blank",
                                 str_detect(pattern = "QC", SampleID) ~ "QC",
                                 TRUE ~ "Sample"))

data_adrc_TIC %>%
  dplyr::mutate(Order = seq_len(n())) %>% # fake order cause it was not provided
  ggscatter("Order", "TIC", add = "reg.line", color = "Type") +
  stat_cor() 

data_t1000_TIC %>%
  dplyr::mutate(Order = seq_len(n())) %>% # fake order cause it was not provided
  ggscatter("Order", "TIC", add = "reg.line", color = "Type") +
  stat_cor()

data_mars_TIC %>%
  dplyr::mutate(Order = seq_len(n())) %>% # fake order cause it was not provided
  ggscatter("Order", "TIC", add = "reg.line", color = "Type") +
  stat_cor() 


# Check sample type
sample_tic_adrc <- data_adrc_TIC %>% dplyr::filter(Type == "Sample") %>% summarise(median(TIC))
blank_tic_adrc <- data_adrc_TIC  %>% dplyr::filter(Type == "Blank") %>% summarise(median(TIC))

sample_tic_t1000 <- data_t1000_TIC %>% dplyr::filter(Type == "Sample") %>% summarise(median(TIC))
blank_tic_t1000 <- data_t1000_TIC  %>% dplyr::filter(Type == "Blank") %>% summarise(median(TIC))

sample_tic_mars <- data_mars_TIC %>% dplyr::filter(Type == "Sample") %>% summarise(median(TIC))
blank_tic_mars <- data_mars_TIC  %>% dplyr::filter(Type == "Blank") %>% summarise(median(TIC))


# Remove sample with low or high TIC
sample_adrc <- data_adrc_TIC %>% dplyr::filter(Type == "Sample") %>% 
  dplyr::filter(TIC > blank_tic_adrc$`median(TIC)`)

# Calculate the interquartile range (IQR)
Q1 <- quantile(sample_adrc$TIC, 0.25)
Q3 <- quantile(sample_adrc$TIC, 0.75)
IQR <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Filter the dataframe to remove outliers
sample_adrc_filter <- sample_adrc %>%
  dplyr::filter(TIC >= lower_bound & TIC <= upper_bound)

# Remove sample with low or high TIC
sample_t1000 <- data_t1000_TIC %>% dplyr::filter(Type == "Sample") %>% 
  dplyr::filter(TIC > blank_tic_t1000$`median(TIC)`)

# Calculate the interquartile range (IQR)
Q1 <- quantile(sample_t1000$TIC, 0.25)
Q3 <- quantile(sample_t1000$TIC, 0.75)
IQR <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Filter the dataframe to remove outliers
sample_t1000_filter <- sample_t1000 %>%
  dplyr::filter(TIC >= lower_bound & TIC <= upper_bound)

# Remove sample with low or high TIC
sample_mars <- data_mars_TIC %>% dplyr::filter(Type == "Sample") %>% 
  dplyr::filter(TIC > blank_tic_mars$`median(TIC)`)

# Calculate the interquartile range (IQR)
Q1 <- quantile(sample_mars$TIC, 0.25)
Q3 <- quantile(sample_mars$TIC, 0.75)
IQR <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Filter the dataframe to remove outliers
sample_mars_filter <- sample_mars %>%
  dplyr::filter(TIC >= lower_bound & TIC <= upper_bound)


sample_adrc_filter %>%
  dplyr::mutate(Order = seq_len(n())) %>% # fake order cause it was not provided
  ggscatter("Order", "TIC", add = "reg.line", color = "Type") +
  stat_cor() 

sample_t1000_filter %>%
  dplyr::mutate(Order = seq_len(n())) %>% # fake order cause it was not provided
  ggscatter("Order", "TIC", add = "reg.line", color = "Type") +
  stat_cor()

sample_mars_filter %>%
  dplyr::mutate(Order = seq_len(n())) %>% # fake order cause it was not provided
  ggscatter("Order", "TIC", add = "reg.line", color = "Type") +
  stat_cor() 


# Filter data
data_t1000_filter <- data_t1000 %>% dplyr::filter(SampleID %in% sample_t1000_filter$SampleID)
data_adrc_filter <- data_adrc %>% dplyr::filter(SampleID %in% sample_adrc_filter$SampleID)
data_mars_filter <- data_mars %>% dplyr::filter(SampleID %in% sample_mars_filter$SampleID)


# Extract raw peak
t1000_raw <- data_t1000_filter %>% dplyr::select(SampleID, `3031`)
colnames(t1000_raw)[2] <- "Molecule"
adrc_raw <- data_adrc_filter %>% dplyr::select(SampleID, `7132`)
colnames(adrc_raw)[2] <- "Molecule"
mars_raw <- data_mars_filter %>% dplyr::select(SampleID, `2458`)
colnames(mars_raw)[2] <- "Molecule"

t1000_z <- t1000_raw %>%
  dplyr::filter(Molecule > 0) %>%
  dplyr::mutate(Molecule_z = (Molecule - mean(Molecule, na.rm = TRUE)) / sd(Molecule, na.rm = TRUE))
adrc_z <- adrc_raw %>%
  dplyr::filter(Molecule > 0) %>%
  dplyr::mutate(Molecule_z = (Molecule - mean(Molecule, na.rm = TRUE)) / sd(Molecule, na.rm = TRUE))
mars_z <- mars_raw %>%
  dplyr::filter(Molecule > 0) %>%
  dplyr::mutate(Molecule_z = (Molecule - mean(Molecule, na.rm = TRUE)) / sd(Molecule, na.rm = TRUE))


# Combine raw, relative, z-scores
t1000_mol <- t1000_raw %>%
  left_join(data_t1000_TIC) %>% dplyr::select(-Type) %>%
  dplyr::mutate(Molecule_RA = Molecule/TIC) %>%
  dplyr::mutate(LogMol = log(Molecule + 1)) %>%
  left_join(t1000_z %>% dplyr::select(-Molecule))

adrc_mol <- adrc_raw %>%
  left_join(data_adrc_TIC) %>% dplyr::select(-Type) %>%
  dplyr::mutate(Molecule_RA = Molecule/TIC) %>%
  dplyr::mutate(LogMol = log(Molecule + 1)) %>%
  left_join(adrc_z %>% dplyr::select(-Molecule))

mars_mol <- mars_raw %>%
  left_join(data_mars_TIC) %>% dplyr::select(-Type) %>%
  dplyr::mutate(Molecule_RA = Molecule/TIC) %>%
  dplyr::mutate(LogMol = log(Molecule + 1)) %>%
  left_join(mars_z %>% dplyr::select(-Molecule))


########
# ADRC #
########
adrc_info <- adrc_mol %>% left_join(adrc_meta, by = c("SampleID" = "UCSD"))
colnames(adrc_info) <- gsub("[- ]", "_", colnames(adrc_info))

# I have 617 samples but for 17 I have no metadata
adrc_info_filter <- adrc_info %>% dplyr::filter(!(is.na(NACC.ID))) # 600 samples
adrc_info_filter %>% distinct(NACC.ID) # 493 unique IDs --> some repeated measures

# Keep one sample per individual with oldest age and highest NACCUDSD
adrc_final <- adrc_info_filter %>%
  group_by(NACC.ID) %>%
  filter(NACCUDSD == max(NACCUDSD)) %>%
  filter(NACCAGE == max(NACCAGE)) %>%
  filter(LogMol == max(LogMol)) %>%
  ungroup() %>% distinct(NACC.ID, .keep_all = TRUE) %>%
  dplyr::filter(NACCAGE != 25) %>% # remove one subject who is just 25 yo
  dplyr::filter(`Kit.Number` != 436182) # reported to be excluded

# Variables of interest: SEX, NACCBMI, NACCAGE, RACE
# DIABET - is important! Look only at T2D
# Check HYPERCHO, B12DEF, ARTHTYPE, APNEA, DEP2YRS, ANXIETY

# NACCMOCA - MoCA Total Score - corrected for education
# CRAFTVRS, CRAFTURS, CRAFTDVR, CRAFTDRE, CRAFTDTI
# NACCUDSD - Cognitive status at UDS visit
# NACCALZD - Presumptive etiologic diagnosis of the cognitive disorder - Alzheimer's disease

# Sex
adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggboxplot(x = "SEX", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

# BMI
adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(NACCBMI > 0 & NACCBMI < 200) %>%
  ggscatter(x = "NACCBMI", y = "LogMol", add = "reg.line") + stat_cor() # effect

# Age
adrc_age_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(NACCBMI > 0 & NACCBMI < 200) %>%
  ggscatter(x = "NACCAGE", y = "LogMol", add = "reg.line", alpha = 0.2,
            xlab = "Age (years)", ylab = "Log(Peak Area)", color = "#006ba6",
            title = "ADRC - 490 Individuals") + stat_cor() + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(NACCBMI > 0 & NACCBMI < 200) %>%
  lm(formula = LogMol ~ NACCAGE + NACCBMI + SEX)
summary(model)

# 3 samples with LogMol < 7.1 are removed cause they constantly classify as outliers by Grubbs test

# Diabetes
adrc_diab_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(DIABET %in% c(0,2)) %>% # keep only T2D
  ggboxplot(x = "DIABET", y = "LogMol", add = "jitter",
            add.params = list(color = "DIABET", alpha = 0.3), 
            palette = c("#1b4965", "#62b6cb"), legend = "none",
            xlab = "Type 2 Diabetes", ylab = "Log(Peak Area)",
            title = "T2D") +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(DIABET %in% c(0,2)) %>%
  dplyr::mutate(DIABET = as.numeric(gsub(2,1,DIABET))) %>%
  lm(formula = LogMol ~ DIABET + NACCAGE + SEX)
summary(model)

# Check if association with T2D is independent from Hypercholesterolemia
model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(DIABET %in% c(0,2)) %>%
  dplyr::mutate(DIABET = as.numeric(gsub(2,1,DIABET))) %>%
  dplyr::filter(HYPCHOL %in% c(0,1)) %>% 
  lm(formula = LogMol ~ DIABET + HYPCHOL + NACCAGE + SEX)
summary(model)

# Check MOCA as non continuous (clinician suggestion)
adrc_moca_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCMOCA > 0 & NACCMOCA < 35) %>%
  dplyr::mutate(MOCA_cat = case_when(NACCMOCA > 25 ~ "Normal",
                                     TRUE ~ "Putative MCI")) %>%
  ggboxplot(x = "MOCA_cat", y = "LogMol", add = "jitter",
          add.params = list(color = "MOCA_cat", alpha = 0.3), 
          palette = c("#1b4965", "#62b6cb"), legend = "none",
          xlab = "MoCA - Education Adjusted", ylab = "Log(Peak Area)",
          title = "ADRC - Montreal Cognitive Assessment (MoCA) Score") + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCMOCA > 0 & NACCMOCA < 35) %>%
  dplyr::mutate(MOCA_cat = case_when(NACCMOCA > 25 ~ "Normal",
                                     TRUE ~ "Putative MCI")) %>%
  lm(formula = LogMol ~ MOCA_cat + SEX + NACCAGE + NACCBMI)
summary(model)

# NACCUDSD - Cognitive status 
adrc_cogn_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCUDSD %in% c(1, 3, 4)) %>%
  dplyr::mutate(NACCUDSD_comb = case_when(NACCUDSD == 3 | NACCUDSD == 4 ~ "zImapired",
                                          TRUE ~ "Unimpaired")) %>%
  dplyr::mutate_at("NACCUDSD_comb", as.factor) %>%
  ggboxplot(x = "NACCUDSD_comb", y = "LogMol", add = "jitter",
            add.params = list(color = "NACCUDSD_comb", alpha = 0.3), 
            palette = c("#1b4965", "#62b6cb"), legend = "none",
            ylab = "Log(Peak Area)", title = "Cognitive Status",
            xlab = "Congitive Status") +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_compare_means()

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCUDSD %in% c(1, 3, 4)) %>%
  dplyr::mutate(NACCUDSD_comb = case_when(NACCUDSD == 3 | NACCUDSD == 4 ~ "Imapired",
                                          TRUE ~ "Healthy")) %>%
  dplyr::mutate_at("NACCUDSD_comb", as.factor) %>%
  lm(formula = LogMol ~ NACCUDSD_comb + NACCAGE + NACCBMI + SEX)
summary(model)

# Check NACCALZD
adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCUDSD %in% c(1,2,3,4)) %>%
  group_by(NACCUDSD, NACCALZD) %>%
  summarise(count = n())

# NACCALZD - MCI AD + Dementia AD
adrc_alz_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::mutate(ALZD_cat = case_when(NACCUDSD == 1 & NACCALZD == 8 ~ "Unimpaired",
                                     (NACCUDSD == 3 & NACCALZD == 1) | (NACCUDSD == 4 & NACCALZD == 1) ~ "AD",
                                     TRUE ~ "Other")) %>%
  dplyr::filter(ALZD_cat != "Other") %>%
  ggboxplot(x = "ALZD_cat", y = "LogMol", add = "jitter",
            add.params = list(color = "ALZD_cat", alpha = 0.3), 
            palette = c("#1b4965", "#62b6cb"), legend = "none",
            ylab = "Log(Peak Area)", title = "Presumptive AD",
            xlab = "Presumptive AD") +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_compare_means()

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCBMI > 0 & NACCBMI < 100) %>%
  dplyr::mutate(ALZD_cat = case_when(NACCUDSD == 1 & NACCALZD == 8 ~ "Unimpaired",
                                     (NACCUDSD == 3 & NACCALZD == 1) | (NACCUDSD == 4 & NACCALZD == 1) ~ "AD",
                                     TRUE ~ "Other")) %>%
  dplyr::filter(ALZD_cat != "Other") %>%
  dplyr::mutate(ALZD_cat = factor(ALZD_cat)) %>%
  lm(formula = LogMol ~ NACCALZD + NACCAGE + NACCBMI + SEX)
summary(model)

# CRAFTVRS, CRAFTURS, CRAFTDVR, CRAFTDRE, CRAFTDTI
CRAFTVRS_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTVRS > 0) %>%
  dplyr::filter(CRAFTVRS < 45) %>%
  ggscatter(x = "CRAFTVRS", y = "LogMol", add = "reg.line",
            title = "Craft Story 21 Recall (Immediate), verbatim scoring",
            xlab = "Score", ylab = "Log(Peak Area)",
            alpha = 0.2, color = "#006ba6") + stat_cor() + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

CRAFTURS_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTURS > 0) %>%
  dplyr::filter(CRAFTURS < 30) %>%
  ggscatter(x = "CRAFTURS", y = "LogMol", add = "reg.line",
            title = "Craft Story 21 Recall (Immediate), paraphrase scoring",
            xlab = "Score", ylab = "Log(Peak Area)",
            alpha = 0.2, color = "#006ba6") + stat_cor() + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

CRAFTDVR_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTDVR > 0) %>%
  dplyr::filter(CRAFTDVR < 45) %>%
  ggscatter(x = "CRAFTDVR", y = "LogMol", add = "reg.line",
            xlab = "Score", ylab = "Log(Peak Area)", title = "memory",
            alpha = 0.2, color = "#006ba6") + stat_cor() + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

CRAFTDRE_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTDRE > 0) %>%
  dplyr::filter(CRAFTDRE < 30) %>%
  ggscatter(x = "CRAFTDRE", y = "LogMol", add = "reg.line",
            title = "Craft Story 21 Recall (Delayed), paraphrase scoring",
            xlab = "Score", ylab = "Log(Peak Area)",
            alpha = 0.2, color = "#006ba6") + stat_cor() + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model1 <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTVRS > 0 & CRAFTVRS < 45) %>%
  lm(formula = LogMol ~ CRAFTVRS + NACCAGE + SEX + NACCBMI)
summary(model1)

model2 <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTDVR > 0 & CRAFTDVR < 45) %>%
  lm(formula = LogMol ~ CRAFTDVR + NACCAGE + SEX + NACCBMI)
summary(model2)

model3 <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTURS > 0 & CRAFTURS < 30) %>%
  lm(formula = LogMol ~ CRAFTURS + NACCAGE + SEX + NACCBMI)
summary(model3)

model4 <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTDRE > 0 & CRAFTDRE < 30) %>%
  lm(formula = LogMol ~ CRAFTDRE + NACCAGE + SEX + NACCBMI)
summary(model4)


#########
# T1000 #
#########
t1000_info <- t1000_mol %>% 
  left_join(t1000_meta_join, by = c("SampleID" = "Specimen Bar Code")) %>%
  dplyr::filter(!(is.na(Age)))
colnames(t1000_info) <- gsub("[- ]", "_", colnames(t1000_info))

# I have 749 samples with metadata
t1000_info %>% distinct(Patient_ID) # 749 unique IDs --> no repeated measures

# Variable of interest: Age, Ethnicity, Gender, InBody_BMI, InBody_DryLeanMass, and InBody_FatMass

# GroupAssignment --> Healthy Control, Mood/Anxiety
# IPAQ_Category --> Level of activity each week scored categorically (Low/Medium/High)

# Sex
t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(Gender))) %>%
  ggboxplot(x = "Gender", y = "LogMol", add = "jitter") + stat_compare_means()

# BMI
t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(InBody_BMI))) %>%
  ggscatter(x = "InBody_BMI", y = "LogMol", add = "reg.line") + stat_cor()

# Age
t1000_age_plot <- t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  ggscatter(x = "Age", y = "LogMol", add = "reg.line", alpha = 0.2,
            xlab = "Age (years)", ylab = "Log(Peak Area)", color = "#43aa8b",
            title = "T1000 - 749 Individuals") + stat_cor() + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model <- t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(Gender))) %>%
  lm(formula = LogMol ~ Age + Gender + InBody_BMI)
summary(model)

# GroupAssignment
t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(GroupAssignment))) %>%
  ggboxplot(x = "GroupAssignment", y = "LogMol", add = "jitter") + stat_compare_means()

# IPAQ_Category
plot_t100_activity <- t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(Gender))) %>%
  dplyr::filter(IPAQ_Category %in% c("Inactive", "Minimally Active", "HEPA Active")) %>%
  ggboxplot(x = "IPAQ_Category", y = "LogMol", add = "jitter", title = "T1000",
            xlab = "Activity Level", ylab = "Log(Peak Area)", legend = "none",
            add.params = list(color = "IPAQ_Category", alpha = 0.5),
            palette = c("#5FBFA2", "#2E7B64", "#204E4A")) +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model <- t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(Gender))) %>%
  dplyr::filter(IPAQ_Category %in% c("Inactive", "Minimally Active", "HEPA Active")) %>%
  dplyr::mutate(IPAQ_Category = factor(levels = c("Inactive", "Minimally Active", "HEPA Active"), IPAQ_Category)) %>%
  lm(formula = LogMol ~ IPAQ_Category + Age + Gender + InBody_DryLeanMass)
summary(model)


########
# MARS #
########
mars_info <- mars_mol %>% 
  left_join(mars_meta_join, by = c("SampleID" = "serumSpecimenBarCode")) %>%
  dplyr::filter(!(is.na(mars_age))) %>%
  dplyr::mutate(mars_age = gsub(">90", "90", mars_age)) %>%
  dplyr::mutate_at("mars_age", as.numeric)
colnames(mars_info) <- gsub("[- ]", "_", colnames(mars_info))

# I have 292 samples with metadata
mars_info %>% distinct(serumMARSID) %>% nrow() # 234 unique IDs

mars_final <- mars_info %>%
  dplyr::filter(diagnosis %in% c("Normal", "MCI-AD", "Dementia-AD")) %>% group_by(serumMARSID) %>%
  dplyr::mutate(diagnosis = factor(levels = c("Normal", "MCI-AD", "Dementia-AD"), diagnosis)) %>%
  dplyr::mutate_at("diagnosis", as.numeric) %>%
  dplyr::filter(diagnosis == max(diagnosis)) %>%
  dplyr::filter(mars_age == max(mars_age)) %>%
  ungroup()

mars_final %>% distinct(reggieid.x) %>% nrow()

# Sex
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(sex))) %>%
  ggboxplot(x = "sex", y = "LogMol", add = "jitter") + stat_compare_means()

# BMI
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggscatter(x = "bmi", y = "LogMol", add = "reg.line") + stat_cor()

# Age
mars_age_plot <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggscatter(x = "mars_age", y = "LogMol", add = "reg.line", alpha = 0.2,
            xlab = "Age (years)", ylab = "Log(Peak Area)", color = "#ff6663",
            title = "MARS - 231 Individuals") + stat_cor() + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(sex))) %>%
  lm(formula = LogMol ~ mars_age + sex + bmi)
summary(model)

# Diagnosis
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::mutate(diagnosis_lump = case_when(diagnosis == 1 ~ "Normal", TRUE ~ "Impaired")) %>%
  ggboxplot(x = "diagnosis_lump", y = "LogMol", add = "jitter") + stat_compare_means() # tendency

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::mutate(diagnosis_lump = case_when(diagnosis == 1 ~ "Normal", TRUE ~ "Impaired")) %>%
  lm(formula = LogMol ~ diagnosis_lump + mars_age + bmi + sex)
summary(model)

# CSF Measurements

# a-synuclein
plot_mars_alpha <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(a_synuclein > 0) %>%
  ggscatter(x = "a_synuclein", y = "LogMol", add = "reg.line", alpha = 0.2,
            xlab = "a_synuclein", color = "#ff6663",
            title = "MARS", ylab = FALSE) + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_cor()

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(a_synuclein > 0) %>%
  lm(formula = LogMol ~ a_synuclein + mars_age + bmi + sex)
summary(model)

# Phosphorylated Tau
plot_ptau_mars <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(p_tau > 0) %>%
  ggscatter(x = "p_tau", y = "LogMol", add = "reg.line", alpha = 0.2,
            xlab = "p_tau", color = "#ff6663",
            title = "MARS", ylab = FALSE) + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_cor()

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(p_tau > 0) %>%
  lm(formula = LogMol ~ p_tau + mars_age + bmi + sex)
summary(model)

# AB42/AB40 Ratio
plot_abratio_mars <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::mutate(ratio = a_beta_1_42/a_beta_1_40) %>%
  dplyr::filter(!is.na(ratio)) %>%
  ggscatter(x = "ratio", y = "LogMol", add = "reg.line", alpha = 0.2,
            xlab = "ratio", color = "#ff6663",
            title = "MARS", ylab = FALSE) + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_cor()

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::mutate(ratio = a_beta_1_42/a_beta_1_40) %>%
  dplyr::filter(!is.na(ratio)) %>%
  lm(formula = LogMol ~ ratio + mars_age + bmi + sex)
summary(model)
