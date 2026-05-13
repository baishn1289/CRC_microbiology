library(tidyverse)
library(ggVennDiagram)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(ggplotify)
library(eulerr)

clean_species <- function(df) {
  df %>%
    mutate(Species = str_remove(feature, ".*s__")) %>% 
    select(Species, coef, pval)
}

df_eo <- clean_species(meta_results_EO)
df_lo <- clean_species(meta_results_LO)
eo_species <- unique(df_eo$Species)
lo_species <- unique(df_lo$Species)
n_eo <- length(eo_species)
n_lo <- length(lo_species)
n_overlap <- length(intersect(eo_species, lo_species))

eo_only <- length(setdiff(eo_species, lo_species))
lo_only <- length(setdiff(lo_species, eo_species))

Jaccard_Index <- round(
  n_overlap / length(unique(c(eo_species, lo_species))),
  3
)
fit <- euler(c(
  "EOCRC" = eo_only,
  "LOCRC" = lo_only,
  "EOCRC&LOCRC" = n_overlap
))

venn_title <- paste0(
  "Overlap of significant CRC-associated species\n",
  "Jaccard index = ", Jaccard_Index
)

p_venn_raw <- plot(fit, fills = list( fill = c("#C83E4D", "#3F5F8F"),alpha = 0.80),
  edges = list( col = "grey25", lwd = 1),
  labels = FALSE,  
  quantities = list(
    fontsize = 13,
    fontface = "plain",
    col = "black"
  )
)

p_venn_gg <- ggplotify::as.ggplot(p_venn_raw) +
  labs(title = venn_title) +
  annotate(
    "text",
    x = 0, y = 0.7,
    label = paste0("EOCRC (n = ", n_eo, ")"),
    fontface = "bold",
    size = 4.2,
    hjust = 0
  ) +
  annotate(
    "text",
    x = 0.45, y = 0.7,
    label = paste0("LOCRC (n = ", n_lo, ")"),
    fontface = "bold",
    size = 4.2,
    hjust = 0
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5,
      lineheight = 0.95,
      margin = margin(b = 6)
    ),
    plot.margin = margin(t = 12, r = 8, b = 8, l = 8)
  )
p_venn_gg

pdf(width = 6, height = 7,
    paste0(base_symbol,
           "downstream_analysis/output/picture/5.3.Venn.pdf"))
print(p_venn_gg)
dev.off()



df_compare <- full_join(df_eo, df_lo, by = "Species", suffix = c("_EO", "_LO")) %>%
  mutate(Status = case_when(
    is.na(coef_LO) ~ "EO Unique",
    is.na(coef_EO) ~ "LO Unique",
    sign(coef_EO) == sign(coef_LO) ~ "Shared (Consistent)",
    TRUE ~ "Shared (Inconsistent)"
  ))


p_compare <- ggplot(df_compare %>% filter(!is.na(coef_EO) & !is.na(coef_LO)), 
                    aes(x = coef_EO, y = coef_LO)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(color = Status), size = 3) +
  geom_text_repel(aes(label = Species), size = 3, max.overlaps = 15) +
  scale_color_manual(values = c("Shared (Consistent)" = "#E41A1C", 
                                "Shared (Inconsistent)" = "#377EB8")) +
  labs(title = "Effect Size Consistency for Shared Species",
       x = "Effect Size in EOCRC (Log Odds/Coef)",
       y = "Effect Size in LOCRC (Log Odds/Coef)") +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid = element_blank(),
    legend.position = "none"
    )   

print(p_compare)

ggsave(filename = "5.3.effect_size_compare.pdf",
       plot = p_compare,
       path = paste0(base_symbol,
                     "downstream_analysis/output/picture/"),
       width = 10.5,
       height = 8,
       dpi = 300)


EO_unique_species <- setdiff(df_eo$Species,df_lo$Species)
meta.crc_EO_unique_species <- meta_results_EO %>%  
  mutate(Species = str_remove(str_extract(feature, 's__[A-Za-z0-9_]+'), '^s__')) %>%
  filter(Species %in% EO_unique_species) %>% 
  select(Species,coef)

meta.crc_LO_unique_species <- meta_results_LO %>%  
  mutate(Species = str_remove(str_extract(feature, 's__[A-Za-z0-9_]+'), '^s__')) %>%
  filter(Species %in% EO_unique_species) %>% 
  select(Species,coef)

unique_species_cal <- left_join(meta.crc_EO_unique_species,meta.crc_LO_unique_species,
                                suffix = c('EO','LO') ,by = 'Species' ) %>% 
  mutate(direction = ifelse(coefEO * coefLO <0 ,'diff','same'))


df_shared_test <- df_compare %>%
  filter(!is.na(coef_EO) & !is.na(coef_LO)) %>%
  mutate(
    abs_EO = abs(coef_EO),
    abs_LO = abs(coef_LO),
    Direction = ifelse(coef_EO > 0, "Upregulated (Positive)", "Downregulated (Negative)")
  )

q1_data <- df_shared_test %>% filter(coef_EO > 0)
res_q1 <- wilcox.test(q1_data$abs_EO, q1_data$abs_LO, paired = TRUE)

q3_data <- df_shared_test %>% filter(coef_EO < 0)
res_q3 <- wilcox.test(q3_data$abs_EO, q3_data$abs_LO, paired = TRUE)


stats_summary <- df_shared_test %>%
  group_size()


df_plot <- df_shared_test %>%
  select(Species, Direction, abs_EO, abs_LO) %>%
  pivot_longer(cols = c(abs_EO, abs_LO), names_to = "Group", values_to = "Intensity") %>%
  mutate(Group = factor(Group, levels = c("abs_EO", "abs_LO"), labels = c("EOCRC", "LOCRC")))

my_comparisons <- list(c("EOCRC", "LOCRC"))

p_strength <- ggplot(df_plot, aes(x = Group, y = Intensity, color = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_line(aes(group = Species), color = "grey", alpha = 0.4) + # 连线显示配对趋势
  geom_point(size = 2) +
  facet_wrap(~Direction) +
  # 核心修改部分
  stat_compare_means(
    comparisons = my_comparisons,          
    paired = TRUE, 
    label = "p.format",                      
    method = "wilcox.test",
    method.args = list(exact = TRUE, correct = FALSE),
    tip.length = 0.02,                      
    label.y = max(df_plot$Intensity) * 1.05   
  ) + 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),    
    plot.subtitle = element_text(hjust = 0.5, size = 12),             
    strip.text = element_text(size = 12, face = "bold"),               
    axis.title.y = element_text(size = 12),                      
    axis.text = element_text(size = 11),                               
    legend.position = "none",                              
    strip.background = element_blank(),
    panel.grid = element_blank()
    ) +
  labs(title = "Comparison of Effect Size Magnitude",
       subtitle = "Paired comparison (Exact Wilcoxon) for shared species",
       y = "Effect Size Magnitude (|Log Odds|)",
       x = "") +
  scale_color_manual(values = c("EOCRC" = "#E41A1C", "LOCRC" = "#377EB8")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))


print(p_strength)
ggsave(paste0(base_symbol,
              "downstream_analysis/output/picture/5.3.Effect_Size_Magnitude_Comparison.pdf"), 
       p_strength, width = 8, height = 6)


layout_design <- c(
  area(t = 1, l = 1, b = 1, r = 4),  
  area(t = 2, l = 1, b = 2, r = 4), 
  area(t = 1, l = 5, b = 2, r = 9)  
)

main_figure2_composite <- p_venn_gg + 
  p_strength + 
  p_compare +
  plot_layout(
    design = layout_design,
    heights = c(1,1)
  ) +
  plot_annotation(
    tag_levels = list(c("C", "D", "E"))
  ) &
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    plot.tag.position = c(0.02, 1)
  )
main_figure2_composite
ggsave(paste0(base_symbol,
              "downstream_analysis/output/picture/5.3.figure2_composite.pdf"), 
       main_figure2_composite, width = 14, height =7)
