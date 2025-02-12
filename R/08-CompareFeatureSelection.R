
# Plot approaches for feature selection. Shows how EWAS and aligning the genome 
# can impact clock accuracy.

#  1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(cowplot)

# Load outcomes from feature selection tests and number of features
comp_fs <- readRDS('output/comp_fs_temp.rds') %>%
  left_join(readRDS('output/n_sites_fs_temp.rds')) %>%
  # Factor EWAS types
  mutate(type = factor(type, 
                       levels = c('meth_dat_norm', 'meth_dat_align', 
                                  'meth_ewas_sex', 'meth_ewas_tis', 
                                  'meth_ewas_tis_sex', 'meth_dat_full_ewas'),
                       labels = c('No F.S.', 'Align.',
                                  'Age, sex', 'Age, tiss.', 
                                  'Age, tiss., sex', 'Full F.S.'))) %>%
  # Calculate mean accuracy by type
  group_by(type) %>%
  mutate(mean_mae = mean(mae), mean_rsq = mean(rsq))

n_fs <- comp_fs %>%
  pivot_longer(cols = c(mae, rsq), names_to = 'metric', values_to = 'ypos') %>%
  group_by(type, metric) %>%
  summarize(ypos = max(ypos), number_sites = unique(number_sites)) %>%
  mutate(ypos = ifelse(metric == 'mae', ypos + 0.1, ypos + 0.01))

# 2 - Boxplots ====

# MAE
mae_plt <- comp_fs %>%
  ggplot(aes(x = type, y = mae)) +
  geom_boxplot(aes(fill = mean_mae)) +
  geom_text(data = n_fs[n_fs$metric == 'mae' ,], 
            aes(label = number_sites, y = ypos),
            position = position_dodge(width = .75),
            show.legend = FALSE, size = 5) +
  scale_fill_gradient(low = '#3388ff', high = '#002255') +
  ylab('Median absolute error') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 1), 'cm'),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, colour = 'black'),
        axis.text.x = element_text(size = 15, colour = 'black', angle = 90, hjust = 1),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black'))

rsq_plt <- comp_fs %>%
  ggplot(aes(x = type, y = rsq)) +
  geom_boxplot(aes(fill = mean_rsq)) +
  geom_text(data = n_fs[n_fs$metric == 'rsq' ,], 
            aes(label = number_sites, y = ypos),
            position = position_dodge(width = 0.75),
            show.legend = FALSE, size = 5) +
  scale_fill_gradient(low = '#552200', high = '#ff5e00') +
  ylab('R-squared') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 1), 'cm'),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, colour = 'black'),
        axis.text.x = element_text(size = 15, colour = 'black', angle = 90, hjust = 1),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black'))

# Make panel plot
plot_grid(rsq_plt, mae_plt, 
          ncol = 2, labels = c('A', 'B'), align = 'v', label_size = 20)

# 3 - Save figure ====

ggsave('figures/pb_data_feature_selection.tiff', plot = last_plot(), 
       device = 'tiff', dpi = 300, height = 12, width = 27, units = 'cm')

