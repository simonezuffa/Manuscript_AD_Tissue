setwd("~/OneDrive - University of California, San Diego Health/Projects/Caltech/Manuscript_AD_Tissue")

library(tidyverse)
library(ggpubr)

metadata_3xtg <- read_csv("data/microbiome/table_plot/3xtg/meta_3xtg_spf.csv")
metadata_5xfad <- read_csv("data/microbiome/table_plot/5xfad/meta_5xfad_spf.csv")

##############
# phylo-RPCA #
##############

phylorpca_3xtg <- read.delim("data/microbiome/table_plot/3xtg/phyloRPCA_ordination_reformat.txt") %>%
  left_join(metadata_3xtg, by = c("X" = "sample_name"))

i <- "study_type"

phylorpca_3xtg_plot <- phylorpca_3xtg %>%
  ggscatter(x = "Axis1_62.5", y = "Axis2_34.4", color = i, alpha = 0.6,
            title = "phylo-RPCA 3xTg Study", palette = c("#B7E6A5", "#007188"),
            xlab = "PC1 (62.5%)", ylab = "PC2 (34.4%)", 
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = phylorpca_3xtg %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("Axis")), mean), aes(Axis1_62.5, Axis2_34.4, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylorpca_3xtg_plot, filename = "phyloRPCA_3xtg.svg", device = "svg", dpi = "retina", height = 2, width = 2)

phylorpca_3xtg_sac <- read.delim("data/microbiome/table_plot/3xtg/sac-phyloRPCA_ordination_reformat.txt") %>%
  left_join(metadata_3xtg, by = c("X" = "sample_name"))

i <- "genotype"

phylorpca_3xtg_sac_plot <- phylorpca_3xtg_sac %>%
  ggscatter(x = "Axis1_63.9", y = "Axis2_30.7", color = i, alpha = 0.6,
            title = "phylo-RPCA 3xTg Sacrifice", palette = c("#5F4B8BFF", "#E69A8DFF"),
            xlab = "PC1 (63.9%)", ylab = "PC2 (30.7%)", 
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = phylorpca_3xtg_sac %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("Axis")), mean), aes(Axis1_63.9, Axis2_30.7, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylorpca_3xtg_sac_plot, filename = "phyloRPCA_3xtg_sacrifice.svg", device = "svg", dpi = "retina", height = 2, width = 2)

i <- "sex"

phylorpca_3xtg_sac_sex_plot <- phylorpca_3xtg_sac %>%
  ggscatter(x = "Axis1_63.9", y = "Axis2_30.7", color = i, alpha = 0.6,
            title = "phylo-RPCA 3xTg Sacrifice",
            xlab = "PC1 (63.9%)", ylab = "PC2 (30.7%)", 
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = phylorpca_3xtg_sac %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("Axis")), mean), aes(Axis1_63.9, Axis2_30.7, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylorpca_3xtg_sac_sex_plot, filename = "phyloRPCA_3xtg_sex_sacrifice.svg", device = "svg", dpi = "retina", height = 2, width = 2)


phylorpca_5xfad <- read.delim("data/microbiome/table_plot/5xfad/phyloRPCA_ordination_reformat.txt") %>%
  left_join(metadata_5xfad)

i <- "intervention"

phylorpca_5xfad_plot <- phylorpca_5xfad %>%
  ggscatter(x = "Axis1_0.790", y = "Axis2_0.193", color = i, alpha = 0.6,
            title = "phylo-RPCA 5xFAD Study", palette = c("#B7E6A5", "#007188"),
            xlab = "PC1 (79.0%)", ylab = "PC2 (19.3%)", 
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = phylorpca_5xfad %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("Axis")), mean), aes(Axis1_0.790, Axis2_0.193, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylorpca_5xfad_plot, filename = "phyloRPCA_5xfad.svg", device = "svg", dpi = "retina", height = 2, width = 2)

phylorpca_5xfad_sac <- read.delim("data/microbiome/table_plot/5xfad/sac_phyloRPCA_ordination_reformat.txt") %>%
  left_join(metadata_5xfad, by = c("X" = "sample_name"))

i <- "mouse_strain"

phylorpca_5xfad_sac_plot <- phylorpca_5xfad_sac %>%
  ggscatter(x = "Axis1_66.8", y = "Axis2_31.8", color = i, alpha = 0.6,
            title = "phylo-RPCA 5xFAD - Sacrifice", palette = c("#5F4B8BFF", "#E69A8DFF"),
            xlab = "PC1 (66.8%)", ylab = "PC2 (31.8%)", 
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = phylorpca_5xfad_sac %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("Axis")), mean), aes(Axis1_66.8, Axis2_31.8, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylorpca_5xfad_sac_plot, filename = "phyloRPCA_5xfad_sacrifice.svg", device = "svg", dpi = "retina", height = 2, width = 2)

i <- "sex"

phylorpca_5xfad_sac_sex_plot <- phylorpca_5xfad_sac %>%
  ggscatter(x = "Axis1_66.8", y = "Axis2_31.8", color = i, alpha = 0.6,
            title = "phylo-RPCA 5xFAD - Sacrifice",
            xlab = "PC1 (66.8%)", ylab = "PC2 (31.8%)", 
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = phylorpca_5xfad_sac %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("Axis")), mean), aes(Axis1_66.8, Axis2_31.8, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylorpca_5xfad_sac_sex_plot, filename = "phyloRPCA_5xfad_sexsacrifice.svg", device = "svg", dpi = "retina", height = 2, width = 2)


###########
# Shannon #
###########

shannon_3xtg <- read_tsv("data/microbiome/table_plot/3xtg/shannon_alpha-diversity.tsv") 
colnames(shannon_3xtg)[1] <- "sample_name"

shannon_3xtg_plot <- shannon_3xtg %>%
  left_join(metadata_3xtg) %>%
  dplyr::filter(study_type == "Sacrifice") %>%
  dplyr::mutate(genotype = factor(genotype, levels = c("B6", "3XTG"))) %>%
  ggboxplot(x = "genotype", y = "shannon_entropy", add = "jitter", 
            add.params = list(color = "genotype"),
            title = "Alpha Diversity - 3xTg Sacrifice",
            palette = c("#5F4B8BFF", "#E69A8DFF"), legend = "none",
            ylab = "Shannon Entropy", xlab = FALSE) + ylim(1, 6) + 
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

shannon_5xfad <- read_tsv("data/microbiome/table_plot/5xfad/shannon_alpha-diversity.tsv") 
colnames(shannon_5xfad)[1] <- "sample_name"

shannon_5xfad_plot <- shannon_5xfad %>%
  left_join(metadata_5xfad) %>%
  dplyr::filter(intervention == "Sacrifice") %>%
  ggboxplot(x = "mouse_strain", y = "shannon_entropy", add = "jitter", 
            add.params = list(color = "mouse_strain"),
            title = "Alpha Diversity - 5xFAD Sacrifice",
            palette = c("#5F4B8BFF", "#E69A8DFF"), legend = "none",
            ylab = "Shannon Entropy", xlab = FALSE) + ylim(1, 6) + 
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

shannon_combined <- ggarrange(shannon_3xtg_plot, shannon_5xfad_plot, nrow = 1)

#ggsave(plot = shannon_combined, filename = "alpha_div.svg", device = "svg", dpi = "retina", height = 2, width = 3)


##################
# Log Ratio 3xTg #
##################

phylo_3xtg_ranks <- read_tsv("data/microbiome/table_plot/3xtg/3XTG_RPCA_topbot4_lr.tsv") 

phylo_3xtg_sac_ranks_plot <- phylo_3xtg_ranks %>%
  left_join(metadata_3xtg, by = c("Sample ID" = "sample_name")) %>%
  dplyr::filter(study_type == "Sacrifice") %>%
  dplyr::mutate(genotype = factor(genotype, levels = c("B6", "3XTG")))%>%
  ggboxplot(x = "genotype", y = "Current_Natural_Log_Ratio", add = "jitter", 
            add.params = list(color = "genotype", alpha = 0.7),
            title = "Log Ratio - 3xTg Sacrifice",
            palette = c("#5F4B8BFF", "#E69A8DFF"), legend = "none",
            ylab = "Ln(Top/Bottom)", xlab = FALSE) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylo_3xtg_sac_ranks_plot, filename = "log_ratio_sac_3xtg.svg", device = "svg", dpi = "retina", height = 2, width = 1.5)

# Longitudinal ratio
phylo_3xtg_long_ranks_plot <- phylo_3xtg_ranks %>%
  left_join(metadata_3xtg, by = c("Sample ID" = "sample_name")) %>%
  dplyr::mutate(genotype = factor(genotype, levels = c("B6", "3XTG")))%>%
  dplyr::mutate(host_age = as.numeric(host_age)) %>%
  dplyr::filter(study_type == "Longitudinal") %>%
  ggscatter(x = "host_age", y = "Current_Natural_Log_Ratio", add = "loess", 
            color = "genotype", title = "Log Ratio - 3xTg Longitudinal", alpha = 0.5,
            palette = c("#5F4B8BFF", "#E69A8DFF"), legend = "none", 
            ylab = "Ln (Top/Bottom)", xlab = "Age (months)") +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylo_3xtg_long_ranks_plot, filename = "log_ratio_sac_long_3xtg.svg", device = "svg", dpi = "retina", height = 2, width = 3)

# Linear mixed effect model
model_3xtg <-  phylo_3xtg_ranks %>%
  left_join(metadata_3xtg, by = c("Sample ID" = "sample_name")) %>%
  dplyr::mutate(genotype = factor(genotype, levels = c("B6", "3XTG")))%>%
  dplyr::mutate(host_age = as.numeric(host_age)) %>%
  dplyr::filter(study_type == "Longitudinal") %>%
  lmerTest::lmer(formula = Current_Natural_Log_Ratio ~ genotype + host_age + (1|host_subject_id))

summary(model_3xtg)

###################
# Log Ratio 5xFAD #
###################

phylo_5xfad_ranks <- read_tsv("data/microbiome/table_plot/5xfad/5XFAD_RPCA_Ax1_topbot10_lr.tsv") 

phylo_5xfad_sac_ranks_plot <- phylo_5xfad_ranks %>%
  left_join(metadata_5xfad, by = c("Sample ID" = "sample_name")) %>%
  dplyr::filter(intervention == "Sacrifice") %>%
  dplyr::mutate(mouse_strain = factor(mouse_strain, levels = c("WT 5XFAD", "Het 5XFAD")))%>%
  ggboxplot(x = "mouse_strain", y = "Current_Natural_Log_Ratio", add = "jitter", 
            add.params = list(color = "mouse_strain", alpha = 0.7),
            title = "Log Ratio - 5xFAD Sacrifice",
            palette = c("#5F4B8BFF", "#E69A8DFF"), legend = "none",
            ylab = "Ln(Top/Bottom)", xlab = FALSE) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylo_5xfad_sac_ranks_plot, filename = "log_ratio_sac_5xfad.svg", device = "svg", dpi = "retina", height = 2, width = 1.5)

# Longitudinal ratio
phylo_5xfad_long_ranks_plot <- phylo_5xfad_ranks %>%
  left_join(metadata_5xfad, by = c("Sample ID" = "sample_name")) %>%
  dplyr::mutate(mouse_strain = factor(mouse_strain, levels = c("WT 5XFAD", "Het 5XFAD")))%>%
  dplyr::mutate(host_age = as.numeric(host_age)) %>%
  dplyr::filter(intervention == "Longitudinal") %>%
  ggscatter(x = "host_age", y = "Current_Natural_Log_Ratio", add = "loess", 
            color = "mouse_strain", title = "Log Ratio - 5xfad Longitudinal", alpha = 0.5,
            palette = c("#5F4B8BFF", "#E69A8DFF"), legend = "none",
            ylab = "Ln (Top/Bottom)", xlab = "Age (months)") +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylo_5xfad_long_ranks_plot, filename = "log_ratio_sac_long_5xfad.svg", device = "svg", dpi = "retina", height = 2, width = 3)

# Linear mixed effect model
model_5xfad <-  phylo_5xfad_ranks %>%
  left_join(metadata_5xfad, by = c("Sample ID" = "sample_name")) %>%
  dplyr::mutate(mouse_strain = factor(mouse_strain, levels = c("WT 5XFAD", "Het 5XFAD")))%>%
  dplyr::mutate(host_age = as.numeric(host_age)) %>%
  dplyr::filter(intervention == "Longitudinal") %>%
  lmerTest::lmer(formula = Current_Natural_Log_Ratio ~ mouse_strain + host_age + (1|host_subject_id))

summary(model_5xfad)

# 5xFAD sex
phylo_5xfad_ranks_sex <- read_tsv("data/microbiome/table_plot/5xfad/5XFAD_phyloRPCA_Ax1_topbot10_lr.tsv") 

phylo_5xfad_sac_sex_ranks_plot <- phylo_5xfad_ranks_sex %>%
  left_join(metadata_5xfad, by = c("Sample ID" = "sample_name")) %>%
  dplyr::filter(intervention == "Sacrifice") %>%
  dplyr::mutate(Current_Natural_Log_Ratio = as.numeric(Current_Natural_Log_Ratio)) %>%
  ggboxplot(x = "sex", y = "Current_Natural_Log_Ratio", add = "jitter", 
            add.params = list(color = "sex", alpha = 0.7),
            title = "Log Ratio - 5xFAD Sacrifice", legend = "none",
            ylab = "Ln(Top/Bottom)", xlab = FALSE) + stat_compare_means() +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = phylo_5xfad_sac_sex_ranks_plot, filename = "log_ratio_sac_5xfad_sex.svg", device = "svg", dpi = "retina", height = 2, width = 1.5)


#############
# ANCOM-BC2 #
#############

# 3xTg
species_3xtg <- read.delim("data/microbiome/table_plot/3xtg/3xtg_species-lvl_table.tsv") %>%
  column_to_rownames("X.OTU.ID") %>% t() %>% as.data.frame()

# Remove species present in less than 5% of samples
prop_zeros <- colMeans(species_3xtg == 0) %>% as.data.frame()
colnames(prop_zeros)[1] <- "Proportion"
cols_to_remove <- prop_zeros %>% dplyr::filter(Proportion > 0.95) %>% rownames_to_column("ID")

species_3xtg_filter <- species_3xtg %>% dplyr::select(-cols_to_remove$ID)

# Remove species which account for less that 1% of total counts
reads_species <- species_3xtg_filter %>% colSums() %>% as.data.frame()
colnames(reads_species)[1] <- "reads"
reads_species_filter <- reads_species %>%
  dplyr::mutate(proportion = reads / sum(reads)) %>%
  dplyr::filter(proportion > 0.00001) %>% rownames_to_column("ID")

species_3xtg_filter2 <- species_3xtg_filter %>% dplyr::select(reads_species_filter$ID)
rownames(species_3xtg_filter2) <- gsub("X14748.", "", rownames(species_3xtg_filter2))

# ANCOM-BC2 - Taxa 3xTg
library(phyloseq)
library(ANCOMBC)

metadata_3xtg_phylo <- metadata_3xtg %>% dplyr::mutate(sample_name = gsub("14748.", "", sample_name)) %>%
  dplyr::filter(study_type == "Sacrifice") %>% dplyr::select(sample_name, genotype, host_age, sex) %>%
  dplyr::filter(sample_name %in% rownames(species_3xtg_filter2)) %>% arrange(sample_name) %>%
  column_to_rownames("sample_name")

species_3xtg_phylo <- species_3xtg_filter2 %>% rownames_to_column("sample_name") %>% arrange(sample_name) %>%
  column_to_rownames("sample_name")

ps_3xtg <- phyloseq(otu_table(species_3xtg_phylo, taxa_are_rows = FALSE), sample_data(metadata_3xtg_phylo))

ancom_output <- ancombc2(data = ps_3xtg, fix_formula = "genotype + sex + host_age",
                         p_adj_method = "holm", pseudo_sens = TRUE, prv_cut = 0.10,
                         group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                         alpha = 0.05, n_cl = 4, verbose = TRUE, global = TRUE,
                         pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                         iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                         em_control = list(tol = 1e-5, max_iter = 100),
                         mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                         trend_control = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                                                              matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
                                                              matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)),
                                              node = list(2, 2, 1), solver = "ECOS", B = 100))

ancom_3xtg_filter <- ancom_output$res %>%
  dplyr::filter(abs(lfc_genotypeB6) > 2.2) %>%
  dplyr::filter(q_genotypeB6 < 0.05) %>%
  dplyr::filter(diff_sexmale == FALSE) %>%
  arrange(desc(lfc_genotypeB6)) %>%
  dplyr::mutate(taxon = gsub(".*;s__", "", taxon)) %>%
  dplyr::mutate(lfc_genotypeB6 = - lfc_genotypeB6) %>%
  dplyr::mutate(Category = case_when(lfc_genotypeB6 > 0 ~ "Enriched 3xTg",
                                     lfc_genotypeB6 < 0 ~ "Depleted 3xTg")) %>%
  mutate(taxon = factor(taxon, levels = taxon))

ancom_3xtg_species <- ancom_3xtg_filter %>% ggplot(aes(x = lfc_genotypeB6, y = taxon, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(xmin = lfc_genotypeB6 - se_genotypeB6, xmax = lfc_genotypeB6 + se_genotypeB6), width = 0.2) +
  scale_fill_manual(values = c("Enriched 3xTg" = "#E69A8DFF", "Depleted 3xTg" = "#5F4B8BFF")) +
  labs(x = "Log Fold Change", y = NULL, fill = NULL, title = "ANCOM-BC2 - 3xTg Sacrifice") +
  theme_pubr(base_size = 14) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.text.y = element_text(hjust = 1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ggsave(plot = ancom_3xtg_species, filename = "ancom_3xtg_species_sacrifice.svg", device = "svg", dpi = "retina", height = 4.5, width = 3)


# ANCOM-BC2 - KO 3xTg
ko_3xtg <- read.delim("data/microbiome/table_plot/3xtg/3xtg_sac_KO_table.tsv") %>%
  column_to_rownames("X.OTU.ID") %>% t() %>% as.data.frame()

# Remove species present in less than 5% of samples
prop_zeros <- colMeans(ko_3xtg == 0) %>% as.data.frame()
colnames(prop_zeros)[1] <- "Proportion"
cols_to_remove <- prop_zeros %>% dplyr::filter(Proportion > 0.95) %>% rownames_to_column("ID")

ko_3xtg_filter <- ko_3xtg %>% dplyr::select(-cols_to_remove$ID)

# Remove species which account for less that 1% of total counts
reads_ko <- ko_3xtg_filter %>% colSums() %>% as.data.frame()
colnames(reads_ko)[1] <- "reads"
reads_ko_filter <- reads_ko %>%
  dplyr::mutate(proportion = reads / sum(reads)) %>%
  dplyr::filter(proportion > 0.00001) %>% rownames_to_column("ID")

reads_ko_filter2 <- ko_3xtg_filter %>% dplyr::select(reads_ko_filter$ID)
rownames(reads_ko_filter2) <- gsub("X14748.", "", rownames(reads_ko_filter2))

ko_3xtg_phylo <- reads_ko_filter2 %>% rownames_to_column("sample_name") %>% arrange(sample_name) %>%
  column_to_rownames("sample_name")

ps_3xtg_ko <- phyloseq(otu_table(ko_3xtg_phylo, taxa_are_rows = FALSE), sample_data(metadata_3xtg_phylo))

ancom_ko_output <- ancombc2(data = ps_3xtg_ko, fix_formula = "genotype + sex + host_age",
                         p_adj_method = "holm", pseudo_sens = TRUE, prv_cut = 0.10,
                         group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                         alpha = 0.05, n_cl = 4, verbose = TRUE, global = TRUE,
                         pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                         iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                         em_control = list(tol = 1e-5, max_iter = 100),
                         mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                         trend_control = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                                                              matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
                                                              matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)),
                                              node = list(2, 2, 1), solver = "ECOS", B = 100))

ancom_3xtg_ko_filter <- ancom_ko_output$res %>%
  dplyr::filter(q_genotypeB6 < 0.05) %>%
  dplyr::filter(abs(lfc_genotypeB6) > 2) %>%
  arrange(desc(lfc_genotypeB6)) %>%
  dplyr::mutate(lfc_genotypeB6 = - lfc_genotypeB6) %>%
  dplyr::mutate(Category = case_when(lfc_genotypeB6 > 0 ~ "Enriched 3xTg",
                                     lfc_genotypeB6 < 0 ~ "Depleted 3xTg")) %>%
  mutate(taxon = factor(taxon, levels = taxon))


# ANCOM-BC2 - Pathway 3xTg
path_3xtg <- read.delim("data/microbiome/table_plot/3xtg/pathway.tsv") %>%
  column_to_rownames("X.OTU.ID") %>% t() %>% as.data.frame() %>%
  rownames_to_column("SampleID") %>% 
  dplyr::filter(str_detect(pattern = "MOUSE", SampleID)) %>% arrange(SampleID) %>%
  column_to_rownames("SampleID")

rownames(path_3xtg) <- gsub("X14748.", "", rownames(path_3xtg))

metadata_3xtg_phylo <- metadata_3xtg %>% dplyr::filter(anonymized_name %in% rownames(path_3xtg)) %>% 
  dplyr::select(anonymized_name, genotype, host_age, sex) %>% arrange(anonymized_name) %>%
  column_to_rownames("anonymized_name")

path_3xtg_phylo <- path_3xtg %>% rownames_to_column("SampleID") %>% 
  dplyr::filter(SampleID %in% rownames(metadata_3xtg_phylo)) %>% arrange(SampleID) %>%
  column_to_rownames("SampleID")

ps_3xtg_path <- phyloseq(otu_table(path_3xtg_phylo, taxa_are_rows = FALSE), sample_data(metadata_3xtg_phylo))

ancom_path_output <- ancombc2(data = ps_3xtg_path, fix_formula = "genotype + sex + host_age",
                            p_adj_method = "holm", pseudo_sens = TRUE, prv_cut = 0.10,
                            group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                            alpha = 0.05, n_cl = 4, verbose = TRUE, global = TRUE,
                            pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                            iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),
                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                            trend_control = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                                                                 matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
                                                                 matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)),
                                                 node = list(2, 2, 1), solver = "ECOS", B = 100))

ancom_3xtg_path_filter <- ancom_path_output$res %>%
  dplyr::filter(q_genotypeB6 < 0.05) %>%
  arrange(desc(lfc_genotypeB6)) %>%
  dplyr::mutate(lfc_genotypeB6 = - lfc_genotypeB6) %>%
  dplyr::mutate(Category = case_when(lfc_genotypeB6 > 0 ~ "Enriched 3xTg",
                                     lfc_genotypeB6 < 0 ~ "Depleted 3xTg")) %>%
  mutate(taxon = factor(taxon, levels = taxon))

ancom_3xtg_path <- ancom_3xtg_path_filter %>% ggplot(aes(x = lfc_genotypeB6, y = taxon, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(xmin = lfc_genotypeB6 - se_genotypeB6, xmax = lfc_genotypeB6 + se_genotypeB6), width = 0.2) +
  scale_fill_manual(values = c("Enriched 3xTg" = "#E69A8DFF", "Depleted 3xTg" = "#5F4B8BFF")) +
  labs(x = "Log Fold Change", y = NULL, fill = NULL, title = "ANCOM-BC2 - 3xTg Sacrifice") +
  theme_pubr(base_size = 14) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.text.y = element_text(hjust = 1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ggsave(plot = ancom_3xtg_path, filename = "ancom_3xtg_path_sacrifice.svg", device = "svg", dpi = "retina", height = 2, width = 2)


# ANCOM-BC2 - Enzyme 3xTg
ec_3xtg <- read.delim("data/microbiome/table_plot/3xtg/ec.tsv") %>%
  column_to_rownames("X.OTU.ID") %>% t() %>% as.data.frame() %>%
  rownames_to_column("SampleID") %>% 
  dplyr::filter(str_detect(pattern = "MOUSE", SampleID)) %>% arrange(SampleID) %>%
  column_to_rownames("SampleID")

# Remove species present in less than 5% of samples
prop_zeros <- colMeans(ec_3xtg == 0) %>% as.data.frame()
colnames(prop_zeros)[1] <- "Proportion"
cols_to_remove <- prop_zeros %>% dplyr::filter(Proportion > 0.95) %>% rownames_to_column("ID")

ec_3xtg_filter <- ec_3xtg %>% dplyr::select(-cols_to_remove$ID)

# Remove species which account for less that 1% of total counts
reads_ec <- ec_3xtg_filter %>% colSums() %>% as.data.frame()
colnames(reads_ec)[1] <- "reads"
reads_ec_filter <- reads_ec %>%
  dplyr::mutate(proportion = reads / sum(reads)) %>%
  dplyr::filter(proportion > 0.00001) %>% rownames_to_column("ID")

reads_ec_filter2 <- ec_3xtg_filter %>% dplyr::select(reads_ec_filter$ID)
rownames(reads_ec_filter2) <- gsub("X14748.", "", rownames(reads_ec_filter2))

ec_3xtg_phylo <- reads_ec_filter2 %>% rownames_to_column("sample_name") %>% 
  arrange(sample_name) %>% column_to_rownames("sample_name")

metadata_3xtg_phylo <- metadata_3xtg %>% dplyr::filter(anonymized_name %in% rownames(ec_3xtg_phylo)) %>% 
  dplyr::select(anonymized_name, genotype, host_age, sex) %>% arrange(anonymized_name) %>%
  column_to_rownames("anonymized_name")

ps_3xtg_ec <- phyloseq(otu_table(ec_3xtg_phylo, taxa_are_rows = FALSE), sample_data(metadata_3xtg_phylo))

ancom_ec_output <- ancombc2(data = ps_3xtg_ec, fix_formula = "genotype + sex + host_age",
                              p_adj_method = "holm", pseudo_sens = TRUE, prv_cut = 0.10,
                              group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 4, verbose = TRUE, global = TRUE,
                              pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                              iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                              trend_control = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                                                                   matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
                                                                   matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)),
                                                   node = list(2, 2, 1), solver = "ECOS", B = 100))

ancom_3xtg_ec_filter <- ancom_ec_output$res %>%
  dplyr::filter(q_genotypeB6 < 0.05) %>%
  dplyr::filter(abs(lfc_genotypeB6) > 1.5) %>%
  arrange(desc(lfc_genotypeB6)) %>%
  dplyr::mutate(lfc_genotypeB6 = - lfc_genotypeB6) %>%
  dplyr::mutate(Category = case_when(lfc_genotypeB6 > 0 ~ "Enriched 3xTg",
                                     lfc_genotypeB6 < 0 ~ "Depleted 3xTg")) %>%
  mutate(taxon = factor(taxon, levels = taxon))


#########
# 5xFAD #
#########

species_5xfad <- read.delim("data/microbiome/table_plot/5xfad/5xfad_species-lvl_table.tsv") %>%
  column_to_rownames("X.OTU.ID") %>% t() %>% as.data.frame()

# Remove species present in less than 5% of samples
prop_zeros <- colMeans(species_5xfad == 0) %>% as.data.frame()
colnames(prop_zeros)[1] <- "Proportion"
cols_to_remove <- prop_zeros %>% dplyr::filter(Proportion > 0.95) %>% rownames_to_column("ID")

species_5xfad_filter <- species_5xfad %>% dplyr::select(-cols_to_remove$ID)

# Remove species which account for less that 1% of total counts
reads_species <- species_5xfad_filter %>% colSums() %>% as.data.frame()
colnames(reads_species)[1] <- "reads"
reads_species_filter <- reads_species %>%
  dplyr::mutate(proportion = reads / sum(reads)) %>%
  dplyr::filter(proportion > 0.00001) %>% rownames_to_column("ID")

species_5xfad_filter2 <- species_5xfad_filter %>% dplyr::select(reads_species_filter$ID)
rownames(species_5xfad_filter2) <- gsub("X14957.", "", rownames(species_5xfad_filter2))

# ANCOM-BC2 - Taxa 5xfad
metadata_5xfad_phylo <- metadata_5xfad %>% dplyr::mutate(sample_name = gsub("14957.", "", sample_name)) %>%
  dplyr::filter(intervention == "Sacrifice") %>% dplyr::select(sample_name, mouse_strain, host_age, sex) %>%
  dplyr::filter(sample_name %in% rownames(species_5xfad_filter2)) %>% arrange(sample_name) %>%
  column_to_rownames("sample_name")

species_5xfad_phylo <- species_5xfad_filter2 %>% rownames_to_column("sample_name") %>% 
  arrange(sample_name) %>% column_to_rownames("sample_name")

ps_5xfad <- phyloseq(otu_table(species_5xfad_phylo, taxa_are_rows = FALSE), sample_data(metadata_5xfad_phylo))

ancom_5xfad_output <- ancombc2(data = ps_5xfad, fix_formula = "mouse_strain + sex + host_age",
                               p_adj_method = "holm", pseudo_sens = TRUE, prv_cut = 0.10,
                               group = "mouse_strain", struc_zero = TRUE, neg_lb = TRUE,
                               alpha = 0.05, n_cl = 4, verbose = TRUE, global = TRUE,
                               pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                               iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                               trend_control = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                                                                    matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
                                                                    matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)),
                                                    node = list(2, 2, 1), solver = "ECOS", B = 100))

ancom_5xfad_filter <- ancom_5xfad_output$res %>%
  dplyr::filter(`q_mouse_strainWT 5XFAD` < 0.05) %>%
  dplyr::mutate(taxon = gsub(".*;s__", "", taxon)) %>%
  dplyr::mutate(lfc_genotypeB6 = - `lfc_mouse_strainWT 5XFAD`) %>%
  dplyr::mutate(Category = case_when(lfc_genotypeB6 > 0 ~ "Enriched 5xfad",
                                     lfc_genotypeB6 < 0 ~ "Depleted 5xfad")) %>%
  mutate(taxon = factor(taxon, levels = taxon))

ancom_5xfad_species <- ancom_5xfad_filter %>% 
  ggplot(aes(x = lfc_genotypeB6, y = taxon, fill = Category)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_errorbar(aes(xmin = lfc_genotypeB6 - `se_mouse_strainWT 5XFAD`, xmax = lfc_genotypeB6 + `se_mouse_strainWT 5XFAD`), width = 0.2) +
  scale_fill_manual(values = c("Enriched 5xfad" = "#E69A8DFF", "Depleted 5xfad" = "#5F4B8BFF")) +
  labs(x = "Log Fold Change", y = NULL, fill = NULL, title = "ANCOM-BC2 - 5xFAD Sacrifice") +
  theme_pubr(base_size = 14) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.text.y = element_text(hjust = 1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ggsave(plot = ancom_5xfad_species, filename = "ancom_5xfad_species_sacrifice.svg", device = "svg", dpi = "retina", height = 2, width = 3)

ancom_5xfad_sex_output <- ancombc2(data = ps_5xfad, fix_formula = "sex + mouse_strain + host_age",
                               p_adj_method = "holm", pseudo_sens = TRUE, prv_cut = 0.10,
                               group = "sex", struc_zero = TRUE, neg_lb = TRUE,
                               alpha = 0.05, n_cl = 4, verbose = TRUE, global = TRUE,
                               pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                               iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                               trend_control = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                                                                    matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
                                                                    matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)),
                                                    node = list(2, 2, 1), solver = "ECOS", B = 100))

ancom_5xfad_sex_filter <- ancom_5xfad_sex_output$res %>%
  dplyr::filter(q_sexmale < 0.05) %>%
  dplyr::mutate(taxon = gsub(".*;s__", "", taxon)) %>%
  dplyr::mutate(lfc_genotypeB6 = - `lfc_mouse_strainWT 5XFAD`) %>%
  dplyr::mutate(Category = case_when(lfc_genotypeB6 > 0 ~ "Enriched 5xfad",
                                     lfc_genotypeB6 < 0 ~ "Depleted 5xfad")) %>%
  mutate(taxon = factor(taxon, levels = taxon))

ancom_5xfad_species <- ancom_5xfad_filter %>% 
  ggplot(aes(x = lfc_genotypeB6, y = taxon, fill = Category)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_errorbar(aes(xmin = lfc_genotypeB6 - `se_mouse_strainWT 5XFAD`, xmax = lfc_genotypeB6 + `se_mouse_strainWT 5XFAD`), width = 0.2) +
  scale_fill_manual(values = c("Enriched 5xfad" = "#E69A8DFF", "Depleted 5xfad" = "#5F4B8BFF")) +
  labs(x = "Log Fold Change", y = NULL, fill = NULL, title = "ANCOM-BC2 - 5xFAD Sacrifice") +
  theme_pubr(base_size = 14) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.text.y = element_text(hjust = 1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ANCOM-BC2 - KO 5xfad
ko_5xfad <- read.delim("data/microbiome/table_plot/5xfad/5xfad_sac_KO_table.tsv") %>%
  column_to_rownames("X.OTU.ID") %>% t() %>% as.data.frame()

# Remove species present in less than 5% of samples
prop_zeros <- colMeans(ko_5xfad == 0) %>% as.data.frame()
colnames(prop_zeros)[1] <- "Proportion"
cols_to_remove <- prop_zeros %>% dplyr::filter(Proportion > 0.95) %>% rownames_to_column("ID")

ko_5xfad_filter <- ko_5xfad %>% dplyr::select(-cols_to_remove$ID)

# Remove species which account for less that 1% of total counts
reads_ko <- ko_5xfad_filter %>% colSums() %>% as.data.frame()
colnames(reads_ko)[1] <- "reads"
reads_ko_filter <- reads_ko %>%
  dplyr::mutate(proportion = reads / sum(reads)) %>%
  dplyr::filter(proportion > 0.00001) %>% rownames_to_column("ID")

ko_5xfad_filter2 <- ko_5xfad_filter %>% dplyr::select(reads_ko_filter$ID)
rownames(ko_5xfad_filter2) <- gsub("X14957.", "", rownames(ko_5xfad_filter2))

ko_5xfad_phylo <- ko_5xfad_filter2 %>% rownames_to_column("sample_name") %>% 
  arrange(sample_name) %>% column_to_rownames("sample_name")

ps_5xfad_ko <- phyloseq(otu_table(ko_5xfad_phylo, taxa_are_rows = FALSE), sample_data(metadata_5xfad_phylo))

ancom_ko_5xfad_output <- ancombc2(data = ps_5xfad_ko, fix_formula = "mouse_strain + sex + host_age",
                                  p_adj_method = "holm", pseudo_sens = TRUE, prv_cut = 0.10,
                                  group = "mouse_strain", struc_zero = TRUE, neg_lb = TRUE,
                                  alpha = 0.05, n_cl = 4, verbose = TRUE, global = TRUE,
                                  pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                                  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                  em_control = list(tol = 1e-5, max_iter = 100),
                                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                                                                       matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
                                                                       matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)),
                                                       node = list(2, 2, 1), solver = "ECOS", B = 100))

ancom_5xfad_ko_filter <- ancom_ko_5xfad_output$res %>%
  dplyr::filter(`q_mouse_strainWT 5XFAD` < 0.05) %>%
  arrange(desc(`lfc_mouse_strainWT 5XFAD`)) %>%
  dplyr::mutate(lfc_genotypeB6 = - `lfc_mouse_strainWT 5XFAD`) %>%
  dplyr::mutate(Category = case_when(lfc_genotypeB6 > 0 ~ "Enriched 5xfad",
                                     lfc_genotypeB6 < 0 ~ "Depleted 5xfad")) %>%
  mutate(taxon = factor(taxon, levels = taxon))


# ANCOM-BC2 - Pathway 5xfad
path_5xfad <- read.delim("data/microbiome/table_plot/5xfad/pathway.tsv") %>%
  column_to_rownames("X.OTU.ID") %>% t() %>% as.data.frame() %>%
  rownames_to_column("SampleID") %>% 
  dplyr::filter(str_detect(pattern = "MOUSE", SampleID)) %>% arrange(SampleID) %>%
  column_to_rownames("SampleID")

rownames(path_5xfad) <- gsub("X14957.", "", rownames(path_5xfad))

metadata_5xfad_phylo <- metadata_5xfad %>% dplyr::mutate(sample_name = gsub("14957.", "", sample_name)) %>%
  dplyr::filter(intervention == "Sacrifice") %>% dplyr::select(sample_name, anonymized_name, mouse_strain, host_age, sex) %>% 
  dplyr::filter(anonymized_name %in% rownames(path_5xfad)) %>% 
  dplyr::select(sample_name, mouse_strain, host_age, sex) %>% arrange(sample_name) %>%
  column_to_rownames("sample_name")

path_5xfad_phylo <- path_5xfad %>% rownames_to_column("SampleID") %>% 
  dplyr::filter(SampleID %in% rownames(metadata_5xfad_phylo)) %>% arrange(SampleID) %>%
  column_to_rownames("SampleID")

ps_5xfad_path <- phyloseq(otu_table(path_5xfad_phylo, taxa_are_rows = FALSE), sample_data(metadata_5xfad_phylo))

ancom_path_output <- ancombc2(data = ps_5xfad_path, fix_formula = "mouse_strain + sex + host_age",
                              p_adj_method = "holm", pseudo_sens = TRUE, prv_cut = 0.10,
                              group = "mouse_strain", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 4, verbose = TRUE, global = TRUE,
                              pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                              iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                              trend_control = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                                                                   matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
                                                                   matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)),
                                                   node = list(2, 2, 1), solver = "ECOS", B = 100))

ancom_5xfad_path_filter <- ancom_path_output$res %>%
  dplyr::filter(`q_mouse_strainWT 5XFAD` < 0.05) %>%
  arrange(desc(`lfc_mouse_strainWT 5XFAD`)) %>%
  dplyr::mutate(lfc_genotypeB6 = - `lfc_mouse_strainWT 5XFAD`) %>%
  dplyr::mutate(Category = case_when(lfc_genotypeB6 > 0 ~ "Enriched 5xfad",
                                     lfc_genotypeB6 < 0 ~ "Depleted 5xfad")) %>%
  mutate(taxon = factor(taxon, levels = taxon))

ancom_5xfad_path <- ancom_5xfad_path_filter %>% ggplot(aes(x = lfc_genotypeB6, y = taxon, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(xmin = lfc_genotypeB6 - `se_mouse_strainWT 5XFAD`, xmax = lfc_genotypeB6 + `se_mouse_strainWT 5XFAD`), width = 0.2) +
  scale_fill_manual(values = c("Enriched 5xfad" = "#E69A8DFF", "Depleted 5xfad" = "#5F4B8BFF")) +
  labs(x = "Log Fold Change", y = NULL, fill = NULL, title = "ANCOM-BC2 - 5xFAD Sacrifice") +
  theme_pubr(base_size = 14) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.text.y = element_text(hjust = 1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ggsave(plot = ancom_5xfad_path, filename = "ancom_5xfad_path_sacrifice.svg", device = "svg", dpi = "retina", height = 2, width = 2)


# ANCOM-BC2 - Enzyme 5xfad
ec_5xfad <- read.delim("data/microbiome/table_plot/5xfad/ec.tsv") %>%
  column_to_rownames("X.OTU.ID") %>% t() %>% as.data.frame() %>%
  rownames_to_column("SampleID") %>% 
  dplyr::filter(str_detect(pattern = "MOUSE", SampleID)) %>% arrange(SampleID) %>%
  column_to_rownames("SampleID")

# Remove species present in less than 5% of samples
prop_zeros <- colMeans(ec_5xfad == 0) %>% as.data.frame()
colnames(prop_zeros)[1] <- "Proportion"
cols_to_remove <- prop_zeros %>% dplyr::filter(Proportion > 0.95) %>% rownames_to_column("ID")

ec_5xfad_filter <- ec_5xfad %>% dplyr::select(-cols_to_remove$ID)

# Remove species which account for less that 1% of total counts
reads_ec <- ec_5xfad_filter %>% colSums() %>% as.data.frame()
colnames(reads_ec)[1] <- "reads"
reads_ec_filter <- reads_ec %>%
  dplyr::mutate(proportion = reads / sum(reads)) %>%
  dplyr::filter(proportion > 0.00001) %>% rownames_to_column("ID")

reads_ec_filter2 <- ec_5xfad_filter %>% dplyr::select(reads_ec_filter$ID)
rownames(reads_ec_filter2) <- gsub("X14957.", "", rownames(reads_ec_filter2))

ec_5xfad_phylo <- reads_ec_filter2 %>% rownames_to_column("sample_name") %>% 
  arrange(sample_name) %>% column_to_rownames("sample_name")

metadata_5xfad_phylo <- metadata_5xfad %>% dplyr::mutate(sample_name = gsub("14957.", "", sample_name)) %>%
  dplyr::filter(intervention == "Sacrifice") %>% dplyr::select(sample_name, anonymized_name, mouse_strain, host_age, sex) %>% 
  dplyr::filter(anonymized_name %in% rownames(ec_5xfad_phylo)) %>% 
  dplyr::select(sample_name, mouse_strain, host_age, sex) %>% arrange(sample_name) %>%
  column_to_rownames("sample_name")

ps_5xfad_ec <- phyloseq(otu_table(ec_5xfad_phylo, taxa_are_rows = FALSE), sample_data(metadata_5xfad_phylo))

ancom_ec_output <- ancombc2(data = ps_5xfad_ec, fix_formula = "mouse_strain + sex + host_age",
                            p_adj_method = "holm", pseudo_sens = TRUE, prv_cut = 0.10,
                            group = "mouse_strain", struc_zero = TRUE, neg_lb = TRUE,
                            alpha = 0.05, n_cl = 4, verbose = TRUE, global = TRUE,
                            pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                            iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),
                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                            trend_control = list(contrast = list(matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE),
                                                                 matrix(c(-1, 0, 1, -1), nrow = 2, byrow = TRUE),
                                                                 matrix(c(1, 0, 1, -1), nrow = 2, byrow = TRUE)),
                                                 node = list(2, 2, 1), solver = "ECOS", B = 100))

ancom_5xfad_ec_filter <- ancom_ec_output$res %>%
  dplyr::filter(`q_mouse_strainWT 5XFAD` < 0.05) %>%
  arrange(desc(`lfc_mouse_strainWT 5XFAD`)) %>%
  dplyr::mutate(lfc_genotypeB6 = - `lfc_mouse_strainWT 5XFAD`) %>%
  dplyr::mutate(Category = case_when(lfc_genotypeB6 > 0 ~ "Enriched 5xfad",
                                     lfc_genotypeB6 < 0 ~ "Depleted 5xfad")) %>%
  mutate(taxon = factor(taxon, levels = taxon))


########
# fldc #
########
fldc <- read_csv("data/microbiome/table_plot/3xtg/3xtg_indv_fldC_summary-counts_sfilt_sac_genotype_by-host_subject_id.csv") %>%
  left_join(metadata_3xtg %>% distinct(host_subject_id, .keep_all = TRUE), by = "host_subject_id")

fldc_plot <- fldc %>%
  dplyr::mutate(`genotype.x` = factor(`genotype.x`, levels = c("B6", "3XTG"))) %>%
  ggboxplot(x = "genotype.x", y = "read_count_matches", add = "jitter", 
            add.params = list(color = "genotype.x", alpha = 0.6),
            title = "bsh - 3xTg Sacrifice", 
            palette = c("#5F4B8BFF", "#E69A8DFF"), legend = "none",
            ylab = "Read Counts", xlab = FALSE) + stat_compare_means() +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8)) +
  coord_flip()

#ggsave(plot = fldc_plot, filename = "fldc_plot.svg", device = "svg", dpi = "retina", width = 3, height = 1)


#######
# bsh #
#######

bsh <- read_csv("data/microbiome/table_plot/3xtg/bsh_summary_count_sac-genotype_tube.csv") %>%
  dplyr::mutate(tube_id2 = as.character(tube_id2)) %>%
  left_join(metadata_3xtg, by = c("tube_id2" = "tube_id"))

bsh_plot <- bsh %>% 
  dplyr::mutate(`genotype.x` = factor(`genotype.x`, levels = c("B6", "3XTG"))) %>%
  ggboxplot(x = "genotype.x", y = "read_count_matches", add = "jitter", 
            add.params = list(color = "genotype.x", alpha = 0.6),
            title = "bsh - 3xTg Sacrifice", 
            palette = c("#5F4B8BFF", "#E69A8DFF"), legend = "none",
            ylab = "Read Counts", xlab = FALSE) + stat_compare_means() +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8)) +
  coord_flip()

#ggsave(plot = bsh_plot, filename = "bsh.svg", device = "svg", dpi = "retina", height = 1.2, width = 2)


##############
# joint-RPCA #
##############

# Extract OGUs from significant species in 3xtg
taxonomy <- read_tsv("data/microbiome/table_plot/3xtg/taxonomy.tsv")

ogu_interest <- taxonomy %>%
  dplyr::filter(str_detect(Taxon, paste(ancom_3xtg_filter$taxon, collapse = "|"))) %>% 
  dplyr::select(1,2)

# joint-RPCA serum
joint_serum <- read_tsv("data/microbiome/jointRPCA/serum.tsv")
serum_metab <- read_csv("data/microbiome/jointRPCA/3xtg_serum_features_info_mn.csv")

serum_metab_interest <- serum_metab %>% 
  dplyr::filter(str_detect(`NPC#class`, patter = "carnitine") | 
                  str_detect(Compound_Name, pattern = regex("carnitine", ignore_case = TRUE))) %>%
  dplyr::filter(!(Feature_ID %in% c("1335", "2127", "2516", "1218", "8443"))) %>%
  dplyr::mutate(Feature = paste("serum", Feature_ID, sep = "_"))

joint_serum_filter <- joint_serum %>% 
  dplyr::filter(featureid %in% ogu_interest$`Feature ID`) %>%
  column_to_rownames("featureid") %>%
  dplyr::select(serum_metab_interest$Feature)

# Export for cytoscaoe
cor_edges <- joint_serum_filter %>%
  as.data.frame() %>%
  rownames_to_column(var = "Microbe") %>% 
  pivot_longer(cols = -Microbe, names_to = "Metabolite", values_to = "Correlation") %>%
  filter(!is.na(Correlation)) %>% 
  dplyr::filter(Metabolite %in% c("serum_3078", "serum_4173", "serum_2821")) %>% 
  left_join(ogu_interest, by = c("Microbe" = "Feature ID")) %>%
  group_by(Taxon, Metabolite) %>%
  summarise(Mean_Correlation = mean(Correlation, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Taxon = gsub(".*; s__", "", Taxon)) %>%
  dplyr::filter(abs(Mean_Correlation) > 0.3) %>%
  dplyr::mutate(abs_cor = abs(Mean_Correlation)) %>%
  dplyr::mutate(value_cor = case_when(Mean_Correlation > 0 ~ "Positive",
                                      TRUE ~ "Negative"))
#write_csv(cor_edges, "correlation_network.csv")
