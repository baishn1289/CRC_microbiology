library(tidyverse)
library(ggVennDiagram)
library(RColorBrewer)
library(ggpubr)

clean_species <- function(df) {
  df %>%
    mutate(Species = str_remove(feature, ".*s__")) %>% 
    select(Species, coef, pval)
}

df_eo <- clean_species(meta_results_EO)
df_lo <- clean_species(meta_results_LO)

venn_list <- list(
  Early_Onset = df_eo$Species,
  Late_Onset = df_lo$Species
)

Jaccard_Index <- round(length(intersect(df_eo$Species, df_lo$Species)) / 
        length(unique(c(df_eo$Species, df_lo$Species))), 3)

p_venn <- ggVennDiagram(
  venn_list, 
  label_alpha = 0, 
  category.names = c("Early-Onset CRC", "Late-Onset CRC")
) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2171B5") +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 2), "cm"),
    text = element_text(size = 14),              
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 13)
  ) + 
  labs(
    title = "Overlap of Significant CRC-associated Species",
    subtitle = paste0("Jaccard Index: ", Jaccard_Index)
  ) + 
  scale_x_continuous(expand = expansion(mult = .2))

pdf(width = 6, height = 7,
    paste0(base_symbol,
           "downstream_analysis/output/picture/5.3.Venn.pdf"))
print(p_venn)
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14))   

print(p_compare)

ggsave(filename = "5.3.effect_size_compare.pdf",
       plot = p_compare,
       path = paste0(base_symbol,
                     "downstream_analysis/output/picture/"),
       width = 10.5,
       height = 8,
       dpi = 300)


EO_unique_species <- setdiff(df_eo$Species,df_lo$Species)
meta.crc_EO_unique_species <- meta_results_EO_csv %>%  
  # mutate(Species = str_remove(str_extract(feature, 's__[A-Za-z0-9_]+'), '^s__')) %>% 
  filter(Species %in% EO_unique_species) %>% 
  select(Species,coef)

meta.crc_LO_unique_species <- meta_results_LO_csv %>%  
  # mutate(Species = str_remove(str_extract(feature, 's__[A-Za-z0-9_]+'), '^s__')) %>% 
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



#  (Paired Boxplot)
df_plot <- df_shared_test %>%
  select(Species, Direction, abs_EO, abs_LO) %>%
  pivot_longer(cols = c(abs_EO, abs_LO), names_to = "Group", values_to = "Intensity") %>%
  mutate(Group = factor(Group, levels = c("abs_EO", "abs_LO"), labels = c("EOCRC", "LOCRC")))

my_comparisons <- list(c("EOCRC", "LOCRC"))

p_strength <- ggplot(df_plot, aes(x = Group, y = Intensity, color = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_line(aes(group = Species), color = "grey", alpha = 0.4) + 
  geom_point(size = 2) +
  facet_wrap(~Direction) +
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
    legend.position = "none"                                          
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
