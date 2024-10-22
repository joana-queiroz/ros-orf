# Plotting putative orf's on sequence

#library(tidyverse)
library(ggplot2)

# Set directory and open dataframe created in rosalind_orf.ipynb
setwd('C:/Users/joanaq/Documents/learning-bioinformatics/projects/ros-orf/level-2')

orf_df <- read.csv('csv-output/sample1kb_putativeorfs.csv')

# Abbreviate ORF IDs in plot RosalindXXXX_ORF_X to ORF_X
orf_df$ID.abbrev <- sub(".*(orf.*)", "\\1", orf_df$Putative.ORF.ID)
# Regex pattern breakdown:
# Match any sequence of characters up until... (.*())
# "orf" and everything that follows (orf.*)

# Re-order ORFs according to location in DNA sequence
orf_df <- orf_df[order(orf_df$ORF.first.bp), ]

# Create the plot with arrows
orf_plot_1 <- ggplot(orf_df) +
  geom_segment(
    aes(x = Start.bp, xend = Stop.bp, y = ID.abbrev, yend = ID.abbrev, color = Strand),
    arrow = arrow(length = unit(0.07, "inches")), size = 1) +
  scale_color_manual(
    limits = c('Sense','Antisense'), 
    values = c("Sense" = "#3bb2e5", "Antisense" = "#e53f4a")) +
  scale_fill_discrete(breaks=c('Sense', 'Antisense')) +
  scale_y_discrete(limits = rev(unique(orf_df$ID.abbrev))) +
  scale_x_continuous(position = "top", n.breaks = 8) +
  labs(x = "DNA bp", y = element_blank(), color = 'Strand') +
  theme(text = element_text(family="mono"),
    axis.line.x = element_line(color = "#3e3c3c49", linewidth = 2, linetype = 1),
    panel.grid.major.y = element_line(color = "#0101015f", linetype = 3),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "#fdfdfd"),
    legend.position = c(0.8,0.8),
    legend.box.background = element_rect(color='black'))
ggsave(path = 'plots', filename = "orf_plot_1.png", plot = orf_plot_1)

library('dplyr')

orf_df_filter1 <- orf_df %>%
  filter(Length > 24)

  orf_plot_2 <- ggplot(orf_df_filter1) +
  geom_segment(
    aes(x = Start.bp, xend = Stop.bp, y = ID.abbrev, yend = ID.abbrev, color = Strand),
    arrow = arrow(length = unit(0.07, "inches")), size = 1) +
  scale_color_manual(
    limits = c('Sense','Antisense'), 
    values = c("Sense" = "#3bb2e5", "Antisense" = "#e53f4a")) +
  scale_fill_discrete(breaks=c('Sense', 'Antisense')) +
  scale_y_discrete(limits = rev(unique(orf_df_filter1$ID.abbrev))) +
  scale_x_continuous(position = "top", n.breaks = 8) +
  labs(x = "DNA bp", y = element_blank(), color = 'Strand') +
  theme(text = element_text(family="mono"),
    axis.line.x = element_line(color = "#3e3c3c49", linewidth = 2, linetype = 1),
    panel.grid.major.y = element_line(color = "#0101015f", linetype = 3),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "#fdfdfd"),
    legend.position = c(0.8,0.8),
    legend.box.background = element_rect(color='black'))
ggsave(path = 'plots', filename = "orf_plot_2.png", plot = orf_plot_2)