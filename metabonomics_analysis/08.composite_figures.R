
library(cowplot)

# Main Figure B-E layout geometry.
left_margin <- 0.02
right_margin <- 0.02
top_margin <- 0.02
bottom_margin <- 0.02
col_gap <- 0.015
row_gap <- 0.04

# Top row uses B:C:D width ratio 3:4:6.
top_ratios <- c(3, 4, 6)

usable_width_top <- 1 - left_margin - right_margin - 2 * col_gap
unit_top <- usable_width_top / sum(top_ratios)

w_B <- top_ratios[1] * unit_top
w_C <- top_ratios[2] * unit_top
w_D <- top_ratios[3] * unit_top

x_B <- left_margin
x_C <- x_B + w_B + col_gap
x_D <- x_C + w_C + col_gap

# Two-row layout; bottom support strip is taller.
row_ratios <- c(1, 1.8)
usable_height <- 1 - top_margin - bottom_margin - row_gap

h_top <- usable_height * row_ratios[1] / sum(row_ratios)
h_bottom <- usable_height * row_ratios[2] / sum(row_ratios)

y_bottom <- bottom_margin
y_top <- y_bottom + h_bottom + row_gap

# Panel E spans 10/14 of the usable top-row width.
w_E <- usable_width_top * 10 / 14
x_E <- left_margin


# Absolute cowplot coordinates.
top_panel_y <- 0.550
top_panel_label_gap <- 0.035
top_panel_label_y <- top_panel_y + h_top + top_panel_label_gap
e_panel_height_scale <- 0.75
e_panel_height <- h_bottom * e_panel_height_scale
e_panel_y <- 0.040
e_label_gap <- 0.020
e_panel_label_y <- e_panel_y + e_panel_height + e_label_gap
e_panel_x <- x_E + 0.015
main_figure_composite <- ggdraw() +
  draw_plot(strict_panel_score_plot,     x = x_B, y = top_panel_y, width = w_B, height = h_top) +
  draw_plot(main_panel_effect_plot,      x = x_C, y = top_panel_y, width = w_C, height = h_top) +
  draw_plot(main_canonical_effect_plot,  x = x_D, y = top_panel_y, width = w_D, height = h_top) +
  draw_plot(support_strip_plot,          x = e_panel_x, y = e_panel_y, width = w_E, height = e_panel_height) +
  draw_plot_label(
    label = c("B", "C", "D", "E"),
    x = c(x_B, x_C, x_D, x_E),
    y = c(top_panel_label_y, top_panel_label_y, top_panel_label_y, e_panel_label_y),
    hjust = -0.05,
    vjust = 1.1,
    fontface = "bold",
    size = 26
  )

save_gg(
  main_figure_composite,
  file.path(plot_dir, "main_figure_nitrogen_diversion_v29.pdf"),
  width =22,
  height = 14.5
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

layout_design <- c(
  area(t = 1, l = 1,  b = 1, r = 2.5),  
  area(t = 2, l = 1,  b = 2, r = 6) 
)
supplementary_figure_composite <- all_panel_member_heatmap/free(extended_panel_score_plot, side = "l")+
  plot_layout(
    design = layout_design,
    heights = c(2,1)
  )+ plot_annotation(
    tag_levels = list(c("A", "B"))
  ) &
  theme(
    plot.tag = element_text(size = 26, face = "bold"),
    plot.tag.position = c(0.02, 1)
  )


save_gg(
  supplementary_figure_composite,
  file.path(plot_dir, "supplementary_figure_composite_v29.pdf"),
  width =15,
  height = 22
)
