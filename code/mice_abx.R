setwd("~/OneDrive - University of California, San Diego Health/Projects/Caltech/Manuscript_AD_Tissue")

library(tidyverse)
library(ggpubr)
library(vegan)
library(patchwork)
library(rstatix)

# Read data and metadata
data <- read_csv("data/mice_abx/gnps_abx.csv")
meta <- read.delim("data/mice_abx/metadata.txt")
meta <- meta %>% dplyr::mutate(filename = gsub(".mzXML", "", filename))

# Data tables
data_mice <- data %>%
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE)

data_mice$SampleID <- gsub(".mzXML Peak area", "", data_mice$SampleID)


# Check TIC 
data_TIC <- data.frame(TIC = rowSums(data_mice %>% column_to_rownames("SampleID"))) %>%
  rownames_to_column("SampleID") %>% arrange(TIC) %>% 
  dplyr::mutate(Type = case_when(str_detect(pattern = "blank", SampleID) ~ "Blank",
                                 str_detect(pattern = "CQ", SampleID) ~ "QC",
                                 TRUE ~ "Sample"))

data_TIC %>%
  dplyr::mutate(Order = seq_len(n())) %>% # fake order
  ggscatter("Order", "TIC", add = "reg.line", color = "Type")


# Remove sample with low or high TIC
sample_filter <- data_TIC %>% dplyr::filter(Type == "Sample") %>%
  dplyr::filter(TIC > 3e9)

Q1 <- quantile(sample_filter$TIC, 0.25)
Q3 <- quantile(sample_filter$TIC, 0.75)
IQR <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

sample_final <- sample_filter %>%
  dplyr::filter(TIC >= lower_bound & TIC <= upper_bound)


# Check carnitines --> (9657 - m/z 280.1545 and 10843 - m/z 294.1703 and 9300 - m/z 266.139)
carnitines <- data_mice %>% dplyr::filter(SampleID %in% sample_final$SampleID) %>%
  dplyr::select("SampleID", "9657", "10843", "9300")

# plot 266
plot_266 <- carnitines %>% left_join(meta, by = c("SampleID" = "filename")) %>%
  dplyr::mutate(LogMol = log(`9300` + 1)) %>%
  ggboxplot(x = "ATTRIBUTE_Traitement", y = "LogMol", add = "jitter",
            xlab = FALSE, title = "m/z 266.139", ylab = "Log(Peak Area)",
            add.params = list(color = "ATTRIBUTE_Traitement", alpha = 0.5), legend = "none",
            palette = c("#214d65", "#E5BF86", "#287DAB")) + ylim(0,15)  +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# plot 280
plot_280 <- carnitines %>% left_join(meta, by = c("SampleID" = "filename")) %>%
  dplyr::mutate(LogMol = log(`9657` + 1)) %>%
  ggboxplot(x = "ATTRIBUTE_Traitement", y = "LogMol", add = "jitter",
            xlab = FALSE, title = "m/z 280.154", ylab = "Log(Peak Area)",
            add.params = list(color = "ATTRIBUTE_Traitement", alpha = 0.5), legend = "none",
            palette = c("#214d65", "#E5BF86", "#287DAB")) + ylim(0,15)  +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# plot 294
plot_294 <- carnitines %>% left_join(meta, by = c("SampleID" = "filename")) %>%
  dplyr::mutate(LogMol = log(`10843` + 1)) %>%
  ggboxplot(x = "ATTRIBUTE_Traitement", y = "LogMol", add = "jitter",
            xlab = FALSE, title = "m/z 294.170", ylab = "Log(Peak Area)",
            add.params = list(color = "ATTRIBUTE_Traitement", alpha = 0.5), legend = "none",
            palette = c("#214d65", "#E5BF86", "#287DAB")) + ylim(0,15)  +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 6))

# Combined plots
plot_carn <- ggarrange(plot_266, plot_280, plot_294, nrow = 1)

#ggsave(plot = plot_carn, filename = "Carn_ABX.svg", device = "svg", dpi = "retina", height = 2, width = 4)

# Stats
carnitines %>% left_join(meta, by = c("SampleID" = "filename")) %>%
  dplyr::mutate(LogMol = log(`9300` + 1)) %>%
  pairwise_wilcox_test(formula = LogMol ~ ATTRIBUTE_Traitement)

carnitines %>% left_join(meta, by = c("SampleID" = "filename")) %>%
  dplyr::mutate(LogMol = log(`9657` + 1)) %>%
  pairwise_wilcox_test(formula = LogMol ~ ATTRIBUTE_Traitement)

carnitines %>% left_join(meta, by = c("SampleID" = "filename")) %>%
  dplyr::mutate(LogMol = log(`10843` + 1)) %>%
  pairwise_wilcox_test(formula = LogMol ~ ATTRIBUTE_Traitement)
