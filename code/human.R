setwd("~/OneDrive - University of California, San Diego Health/Projects/Caltech/Manuscript_AD_Tissue")

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

# My feature of interest has RT ~ 2.20 min in serum of Caltech
# adrc        7132
# t1000       3031
# win         2458


# Read metadata - unfortunatelly I cannot share but can be requested via https://adknowledgeportal.synapse.org/
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

# Variables of interest
# SEX
# NACCBMI
# NACCAGE 
# RACE - mostly white or black
# NACCAPOE - there is variance here
# DIABET - is important! Look only at T2D
# Check HYPERCHO B12DEF ARTHTYPE APNEA DEP2YRS ANXIETY
# Use of drugs, smoke, and other health conditions can be messy

# NACCMMSE - Total MMSE score (using D-L-R-O-W) is not available
# BOSTON - Boston Naming Test (30) - Total score is not available
# NACCMOCA - MoCA Total Score - corrected for education
# CRAFTVRS, CRAFTURS, CRAFTDVR, CRAFTDRE, CRAFTDTI
# NACCGDS - Total GDS Score
# NACCUDSD - Cognitive status at UDS visit
# NACCALZD - Presumptive etiologic diagnosis of the cognitive disorder - Alzheimer's disease
# AMYLPET AMYLCSF FDGAD HIPPATR TAUPETAD CSFTAU FDGFTLD TPETFTLD UDSBENTD
# HD-X pTau181 HD-X pTau217 HD-X NFL HD-X GFAP LumipulseABeta40 LumipulseABeta42

# Sex
adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggboxplot(x = "SEX", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

# Race
adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(RACE %in% c(1,2)) %>%
  ggboxplot(x = "RACE", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

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

# 3 samples with LogMol < 7.1 are removed cause they constantly classofy as outliers by Grubbs test

# APOE Status
adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCAPOE %in% c(1,2,3,4,5,6)) %>%
  dplyr::mutate(ApoRisk = case_when(NACCAPOE == 1 | NACCAPOE == 3 | NACCAPOE == 6 ~ "Neutral",
                                    TRUE ~ "Risk")) %>%
  ggboxplot(x = "ApoRisk", y = "LogMol", add = "jitter") + stat_compare_means() # tendency

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
  glm(formula = DIABET ~ LogMol + NACCAGE + SEX, family = "binomial")
summary(model)

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(DIABET %in% c(0,2)) %>%
  dplyr::mutate(DIABET = as.numeric(gsub(2,1,DIABET))) %>%
  lm(formula = LogMol ~ DIABET + NACCAGE + SEX)
summary(model)

# Hypercholesterolemia
adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(HYPCHOL %in% c(0,1)) %>% 
  ggboxplot(x = "HYPCHOL", y = "LogMol", add = "jitter") + stat_compare_means() # effect

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCAGE > 30) %>%
  dplyr::filter(HYPCHOL %in% c(0,1)) %>% 
  lm(formula = LogMol ~ HYPCHOL + NACCAGE + SEX)
summary(model) # effect

# Check if association with T2D is independent from Hypercholesterolemia
model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(DIABET %in% c(0,2)) %>%
  dplyr::mutate(DIABET = as.numeric(gsub(2,1,DIABET))) %>%
  dplyr::filter(HYPCHOL %in% c(0,1)) %>% 
  glm(formula = DIABET ~ LogMol + NACCAGE + SEX + HYPCHOL, family = "binomial")
summary(model)

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(DIABET %in% c(0,2)) %>%
  dplyr::mutate(DIABET = as.numeric(gsub(2,1,DIABET))) %>%
  dplyr::filter(HYPCHOL %in% c(0,1)) %>% 
  lm(formula = LogMol ~ DIABET + HYPCHOL + NACCAGE + SEX)
summary(model)

# Sleep apnea
adrc_sleep_plot <- adrc_final %>%
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(SLEEPAP %in% c(0,1)) %>% 
  ggboxplot(x = "SLEEPAP", y = "LogMol", add = "jitter",
            add.params = list(color = "SLEEPAP", alpha = 0.3), 
            palette = c("#1b4965", "#62b6cb"), legend = "none",
            xlab = "Sleep Apnea", ylab = "Log(Peak Area)",
            title = "Sleep Apnea") +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(SLEEPAP %in% c(0,1)) %>% 
  glm(formula = SLEEPAP ~ LogMol + NACCAGE + SEX + NACCBMI, family = "binomial")
summary(model)

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(DIABET %in% c(0,2)) %>%
  dplyr::mutate(DIABET = as.numeric(gsub(2,1,DIABET))) %>%
  dplyr::filter(SLEEPAP %in% c(0,1)) %>% 
  lm(formula = LogMol ~ SLEEPAP + DIABET + NACCAGE + SEX)
summary(model)

# Depression
adrc_dep_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCAGE > 30) %>%
  dplyr::filter(DEP %in% c(0,1)) %>% 
  ggboxplot(x = "DEP", y = "LogMol", add = "jitter",
            add.params = list(color = "DEP", alpha = 0.3), 
            palette = c("#1b4965", "#62b6cb"), legend = "none",
            xlab = "Depression", ylab = "Log(Peak Area)",
            title = "Depression") +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCBMI > 0 & NACCBMI < 200) %>%
  dplyr::filter(DEP %in% c(0,1)) %>% 
  glm(formula = DEP ~ LogMol + NACCAGE + NACCBMI + SEX, family = "binomial")
summary(model)

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCBMI > 0 & NACCBMI < 200) %>%
  dplyr::filter(DEP %in% c(0,1)) %>% 
  lm(formula = LogMol ~ DEP + NACCAGE + NACCBMI + SEX)
summary(model)

# B12 Deficiency
adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(VB12DEF %in% c(0,1)) %>% 
  ggboxplot(x = "VB12DEF", y = "LogMol", add = "jitter") + stat_compare_means() # low n

# Arthritis region affected - spine
adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(ARTH %in% c(0, 1)) %>% 
  ggboxplot(x = "ARTH", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

# Anxiety
adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(ANXIET %in% c(0,1)) %>% 
  ggboxplot(x = "ANXIET", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

# PD
adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(PD %in% c(0,1)) %>% 
  ggboxplot(x = "PD", y = "LogMol", add = "jitter") + stat_compare_means() # no affected subjects

# TBI
adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(TBI %in% c(0,1)) %>% 
  ggboxplot(x = "TBI", y = "LogMol", add = "jitter") + stat_compare_means() # no affected subjects

# Seizures
adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(SEIZURES %in% c(0,1)) %>% 
  ggboxplot(x = "SEIZURES", y = "LogMol", add = "jitter") + stat_compare_means() # no affected subjects

# NACCMOCA - MoCA Total Score - corrected for education
adrc_moca_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCMOCA > 0 & NACCMOCA < 35) %>%
  ggscatter(x = "NACCMOCA", y = "LogMol", add = "loess", alpha = 0.2,
            xlab = "MoCA - Education Adjusted", ylab = "Log(Peak Area)", color = "#006ba6",
            title = "ADRC - Montreal Cognitive Assessment (MoCA) Score") + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

library(mgcv)
gam_model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCMOCA > 0 & NACCMOCA < 35) %>%
  gam(LogMol ~ s(NACCMOCA) + SEX + s(NACCAGE) + s(NACCBMI), data = .)
summary(gam_model)

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
  glm(formula = NACCUDSD_comb ~ LogMol + NACCAGE + SEX + NACCBMI, family = "binomial")
summary(model)

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(NACCUDSD %in% c(1, 3, 4)) %>%
  dplyr::mutate(NACCUDSD_comb = case_when(NACCUDSD == 3 | NACCUDSD == 4 ~ "Imapired",
                                          TRUE ~ "Healthy")) %>%
  dplyr::mutate_at("NACCUDSD_comb", as.factor) %>%
  lm(formula = LogMol ~ NACCUDSD_comb + NACCAGE + NACCBMI + SEX)
summary(model)

# Check NACCUDSD and NACCALZD
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
  dplyr::mutate(ALZD_cat = case_when(NACCUDSD == 1 & NACCALZD == 8 ~ "Unimpaired",
                                     (NACCUDSD == 3 & NACCALZD == 1) | (NACCUDSD == 4 & NACCALZD == 1) ~ "AD",
                                     TRUE ~ "Other")) %>%
  dplyr::filter(ALZD_cat != "Other") %>%
  dplyr::mutate(ALZD_cat = factor(ALZD_cat)) %>%
  glm(formula = ALZD_cat ~ LogMol + NACCAGE + NACCBMI + SEX, family = "binomial")
summary(model)

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

# MINT
MINT_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(MINTTOTS > 0 & MINTTOTS < 33) %>%
  ggscatter(x = "MINTTOTS", y = "LogMol", add = "reg.line",
            title = "MINTTOTS",
            xlab = "Score", ylab = "Log(Peak Area)",
            alpha = 0.2, color = "#006ba6") + stat_cor() + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

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

model <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTVRS > 0 & CRAFTVRS < 45) %>%
  lm(formula = LogMol ~ CRAFTVRS + NACCAGE + SEX + NACCBMI)
summary(model)

model <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTDVR > 0 & CRAFTDVR < 45) %>%
  lm(formula = LogMol ~ CRAFTDVR + NACCAGE + SEX + NACCBMI)
summary(model)

model <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTURS > 0 & CRAFTURS < 30) %>%
  lm(formula = LogMol ~ CRAFTURS + NACCAGE + SEX + NACCBMI)
summary(model)

model <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CRAFTDRE > 0 & CRAFTDRE < 30) %>%
  lm(formula = LogMol ~ CRAFTDRE + NACCAGE + SEX + NACCBMI)
summary(model)

# AMYLPET AMYLCSF FDGAD HIPPATR TAUPETAD CSFTAU TPETFTLD UDSBENTD
adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(AMYLPET %in% c(0, 1)) %>% 
  ggboxplot(x = "AMYLPET", y = "LogMol", add = "jitter") + stat_compare_means()

adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(AMYLCSF %in% c(0, 1)) %>% 
  ggboxplot(x = "AMYLCSF", y = "LogMol", add = "jitter") + stat_compare_means()

adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(FDGAD %in% c(0, 1)) %>% 
  ggboxplot(x = "FDGAD", y = "LogMol", add = "jitter") + stat_compare_means()

adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(HIPPATR %in% c(0, 1)) %>% 
  ggboxplot(x = "HIPPATR", y = "LogMol", add = "jitter") + stat_compare_means()

adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(TAUPETAD %in% c(0, 1)) %>% 
  ggboxplot(x = "TAUPETAD", y = "LogMol", add = "jitter") + stat_compare_means()

adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(CSFTAU %in% c(0, 1)) %>% 
  ggboxplot(x = "CSFTAU", y = "LogMol", add = "jitter") + stat_compare_means()

adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(UDSBENTD > 0) %>%
  dplyr::filter(UDSBENTD < 20) %>%
  ggscatter(x = "UDSBENTD", y = "LogMol", add = "reg.line") + stat_cor()

# HD-X pTau181 HD-X pTau217 HD-X NFL HD-X GFAP LumipulseABeta40 LumipulseABeta42
adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(!(is.na(HD_X_pTau181))) %>%
  dplyr::filter(HD_X_pTau181 < 12) %>%
  ggscatter(x = "HD_X_pTau181", y = "LogMol", add = "reg.line") + stat_cor() # no effect

adrc_final %>% 
  dplyr::filter(!(is.na(HD_X_pTau181))) %>%
  dplyr::mutate(ALZD_cat = case_when(NACCUDSD == 1 & NACCALZD == 8 ~ "Unimpaired",
                                     (NACCUDSD == 3 & NACCALZD == 1) | (NACCUDSD == 4 & NACCALZD == 1) ~ "AD",
                                     TRUE ~ "Other")) %>%
  dplyr::filter(ALZD_cat != "Other") %>%
  ggboxplot(x = "ALZD_cat", y = "HD_X_pTau181", add = "jitter") + stat_compare_means()

adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggscatter(x = "HD_X_pTau217", y = "LogMol", add = "reg.line") + stat_cor() # no measures

adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(HD_X_NFL < 60) %>%
  ggscatter(x = "HD_X_NFL", y = "LogMol", add = "reg.line") + stat_cor() # no effect

adrc_final %>% 
  dplyr::filter(HD_X_NFL < 60) %>%
  dplyr::mutate(ALZD_cat = case_when(NACCUDSD == 1 & NACCALZD == 8 ~ "Unimpaired",
                                     (NACCUDSD == 3 & NACCALZD == 1) | (NACCUDSD == 4 & NACCALZD == 1) ~ "AD",
                                     TRUE ~ "Other")) %>%
  dplyr::filter(ALZD_cat != "Other") %>%
  ggboxplot(x = "ALZD_cat", y = "HD_X_NFL", add = "jitter") + stat_compare_means()

adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(HD_X_GFAP < 700) %>%
  ggscatter(x = "HD_X_GFAP", y = "LogMol", add = "reg.line") + stat_cor() # no effect

adrc_final %>% 
  dplyr::filter(HD_X_GFAP < 700) %>%
  dplyr::mutate(ALZD_cat = case_when(NACCUDSD == 1 & NACCALZD == 8 ~ "Unimpaired",
                                     (NACCUDSD == 3 & NACCALZD == 1) | (NACCUDSD == 4 & NACCALZD == 1) ~ "AD",
                                     TRUE ~ "Other")) %>%
  dplyr::filter(ALZD_cat != "Other") %>%
  ggboxplot(x = "ALZD_cat", y = "HD_X_GFAP", add = "jitter") + stat_compare_means()

lum40_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(Lumipulse_ABeta40 < 3000) %>%
  ggscatter(x = "Lumipulse_ABeta40", y = "LogMol", add = "reg.line",
            title = "Lumipulse_ABeta40",
            xlab = "Value", ylab = "Log(Peak Area)",
            alpha = 0.2, color = "#006ba6") +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_cor()

lum42_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(Lumipulse_ABeta42 < 300) %>%
  ggscatter(x = "Lumipulse_ABeta42", y = "LogMol", add = "reg.line",
            title = "Lumipulse_ABeta42",
            xlab = "Value", ylab = "Log(Peak Area)",
            alpha = 0.2, color = "#006ba6") +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_cor()

ratio_lum_plot <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(Lumipulse_ABeta42 < 300) %>%
  dplyr::filter(Lumipulse_ABeta40 < 3000) %>%
  dplyr::mutate(Ratio_Lumi = Lumipulse_ABeta42/Lumipulse_ABeta40) %>%
  ggscatter(x = "Ratio_Lumi", y = "LogMol", add = "reg.line",
            title = "Ratio_Lumi", xlab = "Value", ylab = "Log(Peak Area)",
            alpha = 0.2, color = "#006ba6") + stat_cor() +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(Lumipulse_ABeta40 < 3000) %>%
  lm(formula = LogMol ~ Lumipulse_ABeta40 + NACCAGE + NACCBMI + SEX)
summary(model)

model <- adrc_final %>% 
  dplyr::filter(LogMol > 7.1) %>%
  dplyr::filter(Lumipulse_ABeta42 < 300) %>%
  lm(formula = Lumipulse_ABeta42 ~ LogMol + NACCAGE + NACCBMI + SEX)
summary(model)

model <- adrc_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(Lumipulse_ABeta42 < 300) %>%
  dplyr::filter(Lumipulse_ABeta40 < 3000) %>%
  dplyr::mutate(Ratio_Lumi = Lumipulse_ABeta42/Lumipulse_ABeta40) %>%
  lm(formula = Ratio_Lumi ~ LogMol + NACCAGE + NACCBMI + SEX)
summary(model)


#########
# T1000 #
#########
t1000_info <- t1000_mol %>% left_join(t1000_meta_join, by = c("SampleID" = "Specimen Bar Code")) %>%
  dplyr::filter(!(is.na(Age)))
colnames(t1000_info) <- gsub("[- ]", "_", colnames(t1000_info))

# I have 749 samples with metadata
t1000_info %>% distinct(Patient_ID) # 749 unique IDs --> no repeated measures

# Variable of interest
# Age
# Ethnicity
# Gender
# InBody_BMI
# InBody_DryLeanMass and InBody_FatMass

# GroupAssignment --> Healthy Control, Mood/Anxiety
# IPAQ_Category --> Level of activity each week scored categorically (Low/Medium/High)
# OASIS --> Overall Anxiety Severity And Impairment Scale (OASIS) self-report.
# PHQ --> The Patient Health Questionnaire (PHQ-9) 9-item self-report screening tool 


# Sex
t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(Gender))) %>%
  ggboxplot(x = "Gender", y = "LogMol", add = "jitter") + stat_compare_means() # possible effect

# Race
t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(Ethnicity %in% c("Not Hispanic or Latino", "Hispanic or Latino")) %>%
  ggboxplot(x = "Ethnicity", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

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

# BMI
t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(InBody_BMI))) %>%
  ggscatter(x = "InBody_BMI", y = "LogMol", add = "reg.line") + stat_cor() # no effect

# InBody_DryLeanMass and InBody_FatMass
t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(InBody_DryLeanMass))) %>%
  ggscatter(x = "InBody_DryLeanMass", y = "LogMol", add = "reg.line") + stat_cor() # no effect

t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(InBody_FatMass))) %>%
  ggscatter(x = "InBody_FatMass", y = "LogMol", add = "reg.line") + stat_cor() # no effect

# GroupAssignment
t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(GroupAssignment))) %>%
  ggboxplot(x = "GroupAssignment", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

# BAS_Reward
t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  ggscatter(x = "BAS_Reward", y = "LogMol", add = "reg.line") + stat_cor() # maybe effect

model <- t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  lm(formula = LogMol ~ BAS_Reward + Age + Gender + InBody_BMI)
summary(model) # no effect after correction

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
summary(model) # LogMol higher in Inactive

t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(Diagnosis_1))) %>%
  ggboxplot(x = "Diagnosis_1", y = "LogMol", add = "jitter") + stat_compare_means() + coord_flip()

summary(factor(t1000_info$Diagnosis_1))

t1000_info %>% 
  dplyr::filter(LogMol > 6) %>%
  dplyr::filter(!(str_detect(pattern = "Remission", Diagnosis_1))) %>%
  dplyr::filter(str_detect(pattern = "Healthy|Depressive", Diagnosis_1)) %>%
  dplyr::filter(!(str_detect(pattern = "Psy", Diagnosis_1))) %>%
  dplyr::mutate(groupped = case_when(str_detect(pattern = "Single", Diagnosis_1) ~ "Depression Single",
                                     str_detect(pattern = "Recurrent", Diagnosis_1) ~ "Depression Recurrent",
                                     str_detect(pattern = "High Risk", Diagnosis_1) ~ "Healthy at High Risk",
                                     TRUE ~ "Healthy")) %>%
  dplyr::mutate(groupped = factor(levels = c("Healthy", "Healthy at High Risk", "Depression Single", "Depression Recurrent"), groupped)) %>%
  ggboxplot(x = "groupped", y = "LogMol", add = "jitter", title = "T1000",
            xlab = "Depression", ylab = "Log(Peak Area)", legend = "none",
            add.params = list(color = "groupped", alpha = 0.5),
            palette = c("#5FBFA2", "#4B9F86", "#2E7B64", "#204E4A")) +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

model <- t1000_info %>% 
  dplyr::filter(LogMol > 6) %>%
  dplyr::filter(!(str_detect(pattern = "Remission", Diagnosis_1))) %>%
  dplyr::filter(str_detect(pattern = "Healthy|Depressive", Diagnosis_1)) %>%
  dplyr::filter(!(str_detect(pattern = "Psy", Diagnosis_1))) %>%
  dplyr::mutate(groupped = case_when(str_detect(pattern = "Single", Diagnosis_1) ~ "Depression Single",
                                     str_detect(pattern = "Recurrent", Diagnosis_1) ~ "Depression Recurrent",
                                     str_detect(pattern = "High Risk", Diagnosis_1) ~ "Healthy at High Risk",
                                     TRUE ~ "Healthy")) %>%
  dplyr::mutate(groupped = factor(groupped, levels = c("Healthy", "Healthy at High Risk", "Depression Single", "Depression Recurrent"))) %>%
  lm(formula = LogMol ~ groupped + Age + Gender + InBody_BMI)
summary(model) # LogMol higher in Inactive

t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(str_detect(pattern = "Depressive", Diagnosis_1))) %>%
  dplyr::mutate(groupped = case_when(str_detect(pattern = "Remission", Diagnosis_1) ~ "Remission",
                                     TRUE ~ "Depressed")) %>%
  ggboxplot(x = "Diagnosis_1", y = "LogMol", add = "jitter") + stat_compare_means() + coord_flip()

t1000_info %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(str_detect(pattern = "Healthy", Diagnosis_1)) %>%
  dplyr::mutate(groupped = case_when(str_detect(pattern = "Remission", Diagnosis_1) ~ "Remission",
                                     TRUE ~ "Depressed")) %>%
  ggboxplot(x = "Diagnosis_1", y = "LogMol", add = "jitter") + stat_compare_means()


########
# MARS #
########
mars_info <- mars_mol %>% left_join(mars_meta_join, by = c("SampleID" = "serumSpecimenBarCode")) %>%
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
  ggboxplot(x = "sex", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

# Race
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggboxplot(x = "race", y = "LogMol", add = "jitter") + stat_compare_means() # no samples

# BMI
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggscatter(x = "bmi", y = "LogMol", add = "reg.line") + stat_cor() # small effect

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

# ApoE
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(apoe %in% c("E3_E3", "E3_E4", "E4_E4", "E2_E3")) %>%
  ggboxplot(x = "apoe", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

# ApoE4
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(apoe4))) %>%
  ggboxplot(x = "apoe4", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

# Diagnosis
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggboxplot(x = "diagnosis", y = "LogMol", add = "jitter") + stat_compare_means() # tendency

mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::mutate(diagnosis_lump = case_when(diagnosis == 1 ~ "Normal",
                                           TRUE ~ "Impaired")) %>%
  ggboxplot(x = "diagnosis_lump", y = "LogMol", add = "jitter") + stat_compare_means() # tendency

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(diagnosis != 2) %>%
  lm(formula = LogMol ~ diagnosis + mars_age + bmi + sex)
summary(model)

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::mutate(diagnosis_lump = case_when(diagnosis == 1 ~ "Normal",
                                           TRUE ~ "Impaired")) %>%
  lm(formula = LogMol ~ diagnosis_lump + mars_age + bmi + sex)
summary(model)

# REY
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::rename(REY = `Rey_Auditory_Verbal_Learning_Test_(RAVLT)___Long_Delay_Free_Recall_(Trial_7)_Raw`) %>%
  dplyr::filter(REY > 0) %>%
  ggscatter(x = "REY", y = "LogMol", add = "reg.line") + stat_cor()

# Amyloid positive
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(amyloid_positive))) %>%
  ggboxplot(x = "amyloid_positive", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

# Tau positive
mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!(is.na(pTau_positive))) %>%
  ggboxplot(x = "pTau_positive", y = "LogMol", add = "jitter") + stat_compare_means() # no effect

# Series of CSF measurements
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

plot_ttau_mars <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(t_tau > 0) %>%
  ggscatter(x = "t_tau", y = "LogMol", add = "reg.line", alpha = 0.2,
            xlab = "t_tau", color = "#ff6663",
            title = "MARS", ylab = FALSE) + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_cor()

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  lm(formula = LogMol ~ t_tau + mars_age + bmi + sex)
summary(model)

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

mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggscatter(x = "a_beta_1_40", y = "LogMol", add = "reg.line") + stat_cor() # no effect

mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggscatter(x = "a_beta_1_42", y = "LogMol", add = "reg.line") + stat_cor() # no effect

plot_abratio_mars <-mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::mutate(ratio = a_beta_1_42/a_beta_1_40) %>%
  ggscatter(x = "ratio", y = "LogMol", add = "reg.line", alpha = 0.2,
            xlab = "ratio", color = "#ff6663",
            title = "MARS", ylab = FALSE) + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_cor()

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::mutate(ratio = a_beta_1_42/a_beta_1_40) %>%
  lm(formula = LogMol ~ ratio + mars_age + bmi + sex)
summary(model)

mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!is.na(gfap)) %>%
  ggscatter(x = "gfap", y = "LogMol", add = "reg.line") + stat_cor() # tendency

mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!is.na(il_6)) %>%
  ggscatter(x = "il_6", y = "LogMol", add = "reg.line") + stat_cor() # tendency

mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!is.na(neurogranin)) %>%
  ggscatter(x = "neurogranin", y = "LogMol", add = "reg.line") + stat_cor() # tendency

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!is.na(neurogranin)) %>%
  lm(formula = LogMol ~ neurogranin + mars_age + bmi + sex)
summary(model)

mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!is.na(nfl)) %>%
  ggscatter(x = "nfl", y = "LogMol", add = "reg.line") + stat_cor() # effect

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!is.na(nfl)) %>%
  lm(formula = LogMol ~ nfl + mars_age + bmi + sex)
summary(model)

mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!is.na(s100)) %>%
  ggscatter(x = "s100", y = "LogMol", add = "reg.line") + stat_cor() # no effect

plot_sterm2 <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!is.na(s_trem2)) %>%
  ggscatter(x = "s_trem2", y = "LogMol", add = "reg.line", alpha = 0.2,
            xlab = "s_trem2", color = "#ff6663",
            title = "MARS", ylab = FALSE) + 
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_cor()

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  dplyr::filter(!is.na(s_trem2)) %>%
  lm(formula = LogMol ~ s_trem2 + mars_age + bmi + sex)
summary(model)

plot_ykl40 <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  ggscatter(x = "ykl40",  y = "LogMol", add = "reg.line", alpha = 0.2,
            xlab = "ykl40", color = "#ff6663",
            title = "MARS", ylab = FALSE) +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) + stat_cor()

model <- mars_final %>% 
  dplyr::filter(LogMol > 0) %>%
  lm(formula = LogMol ~ ykl40 + mars_age + bmi + sex)
summary(model)


#########
# PLOTS #
#########

age_plots <- ggarrange(adrc_age_plot, t1000_age_plot, mars_age_plot, nrow = 1)
adrc_plots <- ggarrange(adrc_diab_plot, adrc_dep_plot, adrc_moca_plot, adrc_cogn_plot, adrc_alz_plot, nrow = 1)
mars_plots <- ggarrange(plot_ptau_mars, plot_abratio_mars, plot_mars_alpha, nrow = 1)

#ggsave(plot = age_plots, filename = "u19_age.svg", device = "svg", height = 2, dpi = "retina", width = 5)
#ggsave(plot = adrc_plots, filename = "adrc_box.svg", device = "svg", dpi = "retina", height = 2, width = 5.5)
#ggsave(plot = CRAFTDVR_plot, filename = "Memory_plot.svg", device = "svg", dpi = "retina", height = 2, width = 2)
#ggsave(plot = mars_plots, filename = "mars_plots.svg", device = "svg", dpi = "retina", height = 2, width = 4)


# Combined ROC
mars_comb <- mars_final %>% 
  dplyr::select(LogMol, Molecule_z, diagnosis, mars_age, sex, bmi) %>%
  dplyr::filter(LogMol > 0) %>%
  dplyr::mutate(ALZD_cat = case_when(diagnosis %in% c(2,3) ~ "AD",
                                     TRUE ~ "Unimpaired")) %>%
  dplyr::select(Molecule_z, mars_age, sex, bmi, ALZD_cat) %>%
  dplyr::mutate(sex = case_when(sex == "male" ~ 1,
                                TRUE ~ 2))

adrc_comb <- adrc_final %>% dplyr::select(LogMol, Molecule_z, NACCUDSD, NACCALZD, NACCAGE, NACCBMI, SEX) %>%
  dplyr::filter(LogMol > 7.1) %>% 
  dplyr::mutate(ALZD_cat = case_when(NACCUDSD == 1 & NACCALZD == 8 ~ "Unimpaired",
                                     (NACCUDSD == 3 & NACCALZD == 1) | (NACCUDSD == 4 & NACCALZD == 1) ~ "AD",
                                     TRUE ~ "Other")) %>%
  dplyr::filter(ALZD_cat != "Other") %>%
  dplyr::select(Molecule_z, NACCAGE, SEX, NACCBMI, ALZD_cat) %>%
  dplyr::filter(NACCBMI > 10 & NACCBMI < 50)

colnames(mars_comb) <- c("Mol_z", "age", "sex", "bmi", "dx")
colnames(adrc_comb) <- c("Mol_z", "age", "sex", "bmi", "dx")

combined_data <- rbind(adrc_comb, mars_comb) %>% dplyr::mutate(dx = factor(dx))

# NACCALZD - MCI AD + Dementia AD
model_car <- combined_data %>%
  glm(formula = dx ~ Mol_z + age + sex + bmi, family = "binomial")

# Predicted probabilities
predicted_probs_carn <- predict(model_car, type = "response")

library(pROC)
roc_carn <- roc(combined_data$dx, predicted_probs_carn, ci = TRUE, boot.n = 100)
plot(roc_carn, print.auc = TRUE)
