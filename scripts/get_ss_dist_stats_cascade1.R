library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)

# see https://github.com/druglogics/bool-param-maps/ for more info
# and specifically the script:
# https://github.com/druglogics/bool-param-maps/blob/main/scripts/count_models_ss.R
models_ss_stats = readRDS(file = url('https://raw.githubusercontent.com/druglogics/bool-param-maps/main/data/models_ss_stats.rds'))

models_ss_stats %>% group_by(ss_num) %>% tally() %>%
  ggplot(aes(x = ss_num, y = n, fill = as.factor(ss_num))) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_y_continuous(labels = scales::label_number_si()) +
  geom_text(aes(label = n), vjust = -0.3) +
  geom_text(aes(label = paste0(100 * round(n/nrow(models_ss_stats), digits = 2), "%")), size = 10, vjust = c(2.5, 2.5, -2)) +
  theme_pubr() +
  theme(axis.text.x = element_text(size = 18), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Stable states distribution for CASCADE 1.0 models", x = "Number of stable states", y = "Number of models")
ggsave(filename = 'img/ss_dist_cascade1.png',  dpi = "print", width = 7, height = 5)
