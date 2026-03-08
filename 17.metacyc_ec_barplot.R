
list_succinate <- c(
  "PWY-5677",  "P42-PWY",
  "1.2.7.3", "6.2.1.5",             
  "1.3.5.4" 
)

list_lps <- c(
  "PWY-7316", "PWY-7688", "PWY-6478", "PWY-6953", 
  "2.6.1.33", "6.3.2.37", "6.1.2.1",'6.3.2.49','2.5.1.30'
)

list_polyamines <- c(
  "PWY-6922", "PWY-5004",           
  "ARG+POLYAMINE-SYN", "POLYAMSYN-PWY", 
  "MET-SAM-PWY", "HOMOSER-METSYN-PWY", 
  "PWY-8187", "AST-PWY","PWY-6328",          
  "3.5.3.1",'3.5.3.23',            
  "1.4.1.11",
  "1.4.1.4"
)

ec_list_fig4 <- c(list_succinate, list_lps, list_polyamines)[str_starts(c(list_succinate, list_lps, list_polyamines),'[0-9]')]
get_ec_name_fig4 <- get_ec_info_tidy(ec_list_fig4)
get_ec_name_fig4_final <- get_ec_name_fig4 %>% 
  select(Pathway  = EC_Number,Preferred_Name) %>% 
  mutate(Preferred_Name = ifelse(Pathway =='1.3.5.4',
                                 'Succinate dehydrogenase, formerly fumarate reductase',
                                 Preferred_Name))

all_data_raw <- rbind(metacyc_plot_data_explortary,ec_plot_data_explortary) %>% 
  select(Pathway,Pathway_Short,Group,qval.fdr,exposure,coef,I2,pass_qc) 

length(unique(c(list_succinate, list_lps, list_polyamines)))

table(unique(c(list_succinate, list_lps, list_polyamines)) %in% all_data_raw$Pathway_Short)

c(list_succinate, list_lps, list_polyamines)[!c(list_succinate, list_lps, list_polyamines) %in% all_data_raw$Pathway_Short]


metacyc_ec_interaction_result <- rbind(interaction_results_final,
                                       ec_interaction_results_final %>% 
                                         rename(Pathway_Short = Pathway))

ec_manual_map <- c(
  "2.5.1.30" = "EC 2.5.1.30\nFarnesyl diphosphate synthase"
)

simplify_ec_label_v3 <- function(x) {
  out <- as.character(x)
  
  out <- str_remove_all(
    out,
    "\\s*\\([^\\)]*(deaminating|decarboxylating|forming|ADP-forming|ATP-forming|pyruvate-forming|adding.*units)[^\\)]*\\)"
  )
  
  out <- str_remove_all(
    out,
    ":(\\s*)NAD\\(P\\)\\+|:(\\s*)NAD\\+|:(\\s*)ferredoxin|:(\\s*)quinone"
  )
  
  out <- str_remove_all(
    out,
    "\\b(alpha|beta|gamma|delta|erythro|D|L|\\(2E,6E\\))\\-?"
  )
  
  out <- str_squish(out)

  for (ec in names(ec_manual_map)) {
    out <- ifelse(
      str_detect(out, fixed(ec)),
      ec_manual_map[[ec]],
      out
    )
  }
  
  out
}


simplify_metacyc_label <- function(x) {
  out <- as.character(x)
  is_meta <- str_detect(out, "^(PWY-|[A-Z0-9\\+\\-]+-PWY|AST-PWY|ARG\\+|MET-)")
  if (any(is_meta, na.rm = TRUE)) {
    tmp <- out[is_meta]
    tmp <- str_remove_all(tmp, "\\b(dTDP|UDP|GDP|CMP)-")
    tmp <- str_remove_all(tmp, "\\b(alpha|beta|delta|D|L|Nδ)-")
    tmp <- str_replace(tmp, "O-antigen building blocks biosynthesis", "O-antigen biosynthesis")
    tmp <- str_remove_all(tmp, "\\s*\\([^\\)]*(E\\. coli|putida|P\\. putida|E\\.coli)[^\\)]*\\)")
    
    tmp <- str_squish(tmp)
    out[is_meta] <- tmp
  }
  out
}


plot_data <- left_join(all_data_raw,metacyc_ec_interaction_result,
                       by = 'Pathway_Short') %>% 
  filter(Pathway_Short %in% c(list_succinate, list_lps, list_polyamines)) %>%
  left_join(get_ec_name_fig4_final, by = 'Pathway') %>% 
  mutate(
    Pathway  = ifelse(Pathway =='1.3.5.4','1.3.5.4 (mapped by HUMAnN3)',Pathway),
    Pathway = ifelse(is.na(Preferred_Name),Pathway,paste(Pathway,Preferred_Name,sep = ': ')),
    Category = case_when(
    Pathway_Short %in% list_succinate ~ "Module 1\nReconfigured Carbon Flux\n(Enzymatic Re-weighting of Succinate)",
    Pathway_Short %in% list_lps       ~ "Module 2\nCell Surface Remodeling\n(Hyper-Virulence & Camouflage)",
    Pathway_Short %in% list_polyamines ~ "Module 3\nNitrogen Toxicity\n(Polyamine Blockade & Putrefaction)"
  )) %>%
  mutate(Category = factor(Category, levels = c(
    "Module 1\nReconfigured Carbon Flux\n(Enzymatic Re-weighting of Succinate)",
    "Module 2\nCell Surface Remodeling\n(Hyper-Virulence & Camouflage)",
    "Module 3\nNitrogen Toxicity\n(Polyamine Blockade & Putrefaction)"
  ))) %>%    group_by(Pathway_Short) %>%
  mutate(
    Marker_Sig = case_when(`qval.fdr` < 0.1 ~ "Significant (q < 0.1)",
                           `qval.fdr` >= 0.1 ~ "Non-significant",
                           is.na(`qval.fdr`) ~"NA"
                           ),
    Marker_qc = ifelse(pass_qc, "Consistent (QC_pass)", NA),
    Marker_P_interaction = if_else(
      Nominal_Significant == "Yes (*)" &
        abs(coef) == max(abs(coef), na.rm = TRUE),
      "Specific (P_Interaction < 0.05)",
      NA_character_
    )
  )  %>%  ungroup() %>%   
  mutate(
    ID_type = if_else(
      str_detect(Pathway_Short, "^[A-Za-z]"),
      "metacyc","ec" )) %>%
  mutate(
    Pathway = Pathway %>%
      gsub("&alpha;", "α", ., fixed = TRUE) %>%
      gsub("&delta;", "δ", ., fixed = TRUE)
  ) %>% 
      mutate(
        Pathway = Pathway %>%
          gsub("^([0-9])", "EC \\1", .)%>%
          simplify_ec_label_v3() %>%
          simplify_metacyc_label() %>%
          sub(":", "\n", .)
      # gsub("\\band\\b", "and\n", ., perl = TRUE) %>%
      # gsub("\\bderived\\b", "derived\n", ., perl = TRUE)
      ) 

min_coef <- min(plot_data$coef, na.rm = TRUE)
data_range <- max(plot_data$coef, na.rm = TRUE) - min_coef


offset_sig <- 0.16  
offset_i2  <- 0.12 
offset_p_interaction  <- 0.08

marker_pos_sig <- min_coef - (data_range * offset_sig)
marker_pos_i2  <- min_coef - (data_range * offset_i2)
marker_pos_p_interaction  <- min_coef - (data_range * offset_p_interaction)

order_levels <- plot_data %>%
  filter(Group == "EO") %>%
  group_by(Category) %>% 
  arrange(
    factor(ID_type, levels = c("ec","metacyc")),
    abs(coef),
    .by_group = TRUE)%>% 
  ungroup() %>% pull(Pathway) 

order_levels

plot_data$Pathway <- factor(plot_data$Pathway, levels = order_levels)
  
p <- ggplot(plot_data, aes(y = Pathway, x = coef, group = Group)) +

  geom_col(aes(fill = Group), position = position_dodge(width = 0.8), width = 0.7) +
  
  geom_point(aes(x = marker_pos_sig, color = Marker_Sig), 
             position = position_dodge(width = 0.8), 
             shape = 15, size = 3, na.rm = TRUE) +
  geom_point(aes(x = marker_pos_p_interaction, color = Marker_P_interaction), 
             position = position_dodge(width = 0.8), 
             shape = 15, size = 3, na.rm = TRUE) +
  geom_point(aes(x = marker_pos_i2, color = Marker_qc), 
             position = position_dodge(width = 0.8), 
             shape = 15, size = 3, na.rm = TRUE) +

  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +

  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +

  scale_color_manual(
    name = "Statistical Metrics", 

    breaks = c(
      "Significant (q < 0.1)",
      "Non-significant",
      "NA",
      "Consistent (QC_pass)",
      "Specific (P_Interaction < 0.05)"
    ),

    values = c(
      "Significant (q < 0.1)"            = "#fc8d59", 
      "Non-significant"        = "#a65628", 
      "NA"  = "#d9d9d9",  
      "Consistent (QC_pass)"              = "#984ea3",  
      "Specific (P_Interaction < 0.05)"   = "#8da0cb"   
    ),

    na.translate = FALSE
  ) +
  
  scale_fill_manual(
    name = "Cohort Group",
    values = c("EO" = "#EE0000", "LO" = "#3B4992")
  ) + 
  
  guides(
    fill = guide_legend(order = 1), 
    color = guide_legend(
      order = 2, 
      title.position = "top",
      nrow = 2,
      byrow = TRUE 
    )
  ) +
  theme_bw() +
  labs(x = "Pooled-analysis Coefficient (Effect Size)", 
       y = NULL, 
       title = "Age-Stratified Pooled-analysis of Microbial Pathways in Colorectal Cancer") +
  
  theme(
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11),
    strip.text = element_text(size = 11, face = "bold", lineheight = 0.8), 
    strip.background = element_rect(fill = "grey95"),

    panel.grid = element_blank(), 

    axis.ticks.y = element_blank(),

    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",  
    legend.margin = margin(t = 10),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold")
  ) #+


final_xlim_min <- marker_pos_sig * 1.05
p <- p + coord_cartesian(xlim = c(final_xlim_min, max(plot_data$coef) * 1.05), clip = "off")


print(p)


ggsave(paste0(base_symbol,"downstream_analysis/output/picture/18.pwy_ec_barplot.pdf"), 
       p, width = 12.72727, height = 18,
       device = cairo_pdf)  

