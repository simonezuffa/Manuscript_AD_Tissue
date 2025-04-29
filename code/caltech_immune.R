###################################
# To be run after caltech_serum.R #
###################################

setwd("~/OneDrive - University of California, San Diego Health/Projects/Caltech/Manuscript_AD_Tissue")

# Check carnitines of interest
# 2821 - 266.137
# 3078 - 280.153
# 4173 - 294.169

car_int <- data_frame(SampleID = c(data_3xtg$SampleID),
                      Carnitine1 = c(data_3xtg$`3078`),
                      Carnitine2 = c(data_3xtg$`4173`),
                      Carnitine3 = c(data_3xtg$`2821`)) %>%
  dplyr::mutate(LogCarn1 = log(Carnitine1 + 1)) %>% 
  dplyr::mutate(LogCarn2 = log(Carnitine2 + 1)) %>%
  dplyr::mutate(LogCarn3 = log(Carnitine3 + 1)) %>%
  left_join(metadata_metabolomics)


# Read immunological data
library(readxl)

immune_il17 <- read_excel("data/mice_caltech/immune/3xtg/072624_3xTg_IL-17A.xlsx")
immune_cells <- read_excel("data/mice_caltech/immune/3xtg/121423_3xTg_T_B_CD4_CD8_Freq_Counts.xlsx")

car_immune <- car_int %>% 
  dplyr::mutate(Mouse = gsub("MOUSE", "", sample_name)) %>%
  dplyr::mutate(Mouse = gsub("\\..*", "", Mouse)) %>%
  dplyr::mutate_at("Mouse", as.numeric)

il17_car <- immune_il17 %>% left_join(car_immune)


# Plots
il17_car %>% dplyr::filter(LogCarn1 > 0) %>%
  dplyr::filter(Tissue == "MLN") %>% # correlates only to MLN
  ggscatter(x = "IL17A", y = "LogCarn1", add = "reg.line", scales = "free_x", conf.int = TRUE) + stat_cor() 

il17_car %>% dplyr::filter(LogCarn2 > 0) %>%
  dplyr::filter(Tissue == "MLN") %>% # correlates only to MLN
  ggscatter(x = "IL17A", y = "LogCarn2", add = "reg.line", scales = "free_x", conf.int = TRUE) + stat_cor() 

# Combined plot
il17_car_long <- il17_car %>% dplyr::filter(Tissue == "MLN") %>%
  dplyr::select(Sex, IL17A, LogCarn1, LogCarn2, LogCarn3) %>%
  dplyr::filter(!(is.na(LogCarn1))) %>% 
  pivot_longer(cols = c("LogCarn1", "LogCarn2", "LogCarn3"), names_to = "Carn", values_to = "Value") %>% 
  dplyr::filter(Value > 0)

il17_car_plot <- il17_car_long %>%
  ggscatter(x = "IL17A", y = "Value", add = "reg.line", xlab = "Frequency of IL-17A+ in CD4+ T cells (%)",
            color = "Carn", , ylab = "Log(Peak Area)", alpha = 0.5, legend = "none",
            title = "Proinflammatory Response in Mesenteric Lymph Nodes",
            palette = c("#2F3D70", "#D04E59", "#FAE093")) + ylim(4,13) +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = il17_car_plot, filename = "IL17A_plot.svg", device = "svg", dpi = "retina", height = 2, width = 2.5)

model1 <- il17_car %>% dplyr::filter(Tissue == "MLN") %>%
  dplyr::select(Sex, IL17A, LogCarn1) %>%
  dplyr::filter(!(is.na(LogCarn1))) %>% 
  dplyr::filter(LogCarn1 > 0) %>%
  lm(formula = LogCarn1 ~ IL17A + Sex)

summary(model1)

model2 <- il17_car %>% dplyr::filter(Tissue == "MLN") %>%
  dplyr::select(Sex, IL17A, LogCarn2) %>%
  dplyr::filter(!(is.na(LogCarn2))) %>% 
  dplyr::filter(LogCarn2 > 0) %>%
  lm(formula = LogCarn2 ~ IL17A + Sex)

summary(model2)

model3 <- il17_car %>% dplyr::filter(Tissue == "MLN") %>%
  dplyr::select(Sex, IL17A, LogCarn3) %>%
  dplyr::filter(!(is.na(LogCarn3))) %>% 
  dplyr::filter(LogCarn3 > 0) %>%
  lm(formula = LogCarn3 ~ IL17A + Sex)

summary(model3)


# Cell counts
cells_car <- immune_cells %>% left_join(car_immune) %>%
  dplyr::mutate(CD3_CD19 = T_cell_count/B_cell_count) %>%
  dplyr::mutate(CD4_CD8 = log(CD4_T_cell_count/CD8_T_cell_count + 1))

# Combined plot
cells_car_long <- cells_car %>% dplyr::filter(Tissue == "MLN") %>%
  dplyr::select(Sex, CD4_CD8, LogCarn1, LogCarn2, LogCarn3) %>%
  dplyr::filter(!(is.na(LogCarn1))) %>% 
  pivot_longer(cols = c("LogCarn1", "LogCarn2", "LogCarn3"), names_to = "Carn", values_to = "Value")

cells_car_plot <- cells_car_long %>%
  ggscatter(x = "CD4_CD8", y = "Value", add = "reg.line", xlab = "Log(CD4+/CD8+)",
            color = "Carn", legend = "none", ylab = "Log(Peak Area)", alpha = 0.5,
            title = "Ratio Helper to Cytotoxic T Cells in MLN",
            palette = c("#2F3D70", "#D04E59", "#FAE093")) + ylim(4,13) +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = cells_car_plot, filename = "CD4CD8_plot.svg", device = "svg", dpi = "retina", height = 2, width = 2.5)

model4 <- cells_car %>% dplyr::filter(Tissue == "MLN") %>%
  dplyr::select(Sex, CD4_CD8, LogCarn1) %>%
  dplyr::filter(!(is.na(LogCarn1))) %>% 
  dplyr::filter(LogCarn1 > 0) %>%
  lm(formula = LogCarn1 ~ CD4_CD8 + Sex)

summary(model4)

model5 <- cells_car %>% dplyr::filter(Tissue == "MLN") %>%
  dplyr::select(Sex, CD4_CD8, LogCarn2) %>%
  dplyr::filter(!(is.na(LogCarn2))) %>% 
  dplyr::filter(LogCarn2 > 0) %>%
  lm(formula = LogCarn2 ~ CD4_CD8 + Sex)

summary(model5)

model6 <- cells_car %>% dplyr::filter(Tissue == "MLN") %>%
  dplyr::select(Sex, CD4_CD8, LogCarn3) %>%
  dplyr::filter(!(is.na(LogCarn3))) %>% 
  dplyr::filter(LogCarn3 > 0) %>%
  lm(formula = LogCarn3 ~ CD4_CD8 + Sex)

summary(model6)

combined_plots <- ggarrange(il17_car_plot, cells_car_plot)

#ggsave(plot = combined_plots, filename = "immune_plot.svg", device = "svg", dpi = "retina", height = 2, width = 4.5)
