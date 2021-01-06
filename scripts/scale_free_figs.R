library(dplyr)
library(latex2exp)
library(ggplot2)

# load the abmlog-related stats
tbl1 = readRDS(file = "data/abmlog_stats_gamma2.5.rds")
tbl2 = readRDS(file = "data/abmlog_stats_gamma2.rds")

sf_topo_stats = dplyr::bind_rows(tbl1, tbl2)

sf_topo_stats %>%
  ggplot(aes(x = gamma, y = lo_nodes_num, fill = as.factor(gamma))) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(2,2.5)) +
  guides(fill = guide_legend(title = latex2exp::TeX("$\\gamma$"))) +
  labs(y = "Number of link operator nodes", x = latex2exp::TeX("Scale-free exponent $\\gamma$"),
    title = "") +
  theme_classic(base_size = 14)
ggsave(filename = "img/lo_eq_density.png", dpi = "print", width = 7, height = 5)

sf_topo_stats %>%
  ggplot(aes(x = gamma, y = max_num_reg, fill = as.factor(gamma))) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(breaks = c(2,2.5)) +
  guides(fill = guide_legend(title = latex2exp::TeX("$\\gamma$"))) +
  labs(y = "Maximum in-degree", x = latex2exp::TeX("Scale-free exponent $\\gamma$"),
    title = "") +
  theme_classic(base_size = 14)
ggsave(filename = "img/max_reg_density.png", dpi = "print", width = 7, height = 5)

sf_topo_stats %>%
  filter(ss_total < 200000) %>% # exclude 2 'outlier' points
  ggplot(aes(y = ss_total, x = gamma, fill = as.factor(gamma))) +
  geom_boxplot(outlier.alpha = 0) +
  scale_y_continuous(labels = scales::label_number_si()) +
  scale_x_continuous(breaks = c(2,2.5)) +
  scale_fill_brewer(palette = "Set1") +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  guides(fill = guide_legend(title = latex2exp::TeX("$\\gamma$"))) +
  labs(x = latex2exp::TeX("Scale-free exponent $\\gamma$"), title = "Number of Stable states per abmlog-ensemble", y ="") + # for each abmlog-ensemble generated from a scale-free topology
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "img/ss_dist.png", dpi = "print", width = 7, height = 5)

# scale-free topologies percentage that every possible AND-NOT/OR-NOT link-operator configuration resulted in a zero stable state model
sf_topo_stats %>%
  group_by(gamma) %>%
  summarise(zero_ss_per = sum(ss_total == 0)/n(), .groups = 'drop') %>%
  ggplot(aes(x = gamma, y = zero_ss_per, fill = as.factor(gamma))) +
  geom_col() +
  geom_text(aes(label = scales::percent(zero_ss_per)), vjust = -0.5, size = 10) +
  scale_x_continuous(breaks = c(2,2.5)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = guide_legend(title = latex2exp::TeX("$\\gamma$"))) +
  labs(x = latex2exp::TeX("Scale-free exponent $\\gamma$"),
    title = "Zero stable states for all abmlog models",
    y = "Percentage of scale-free topologies") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "img/zero_ss_dist.png", dpi = "print", width = 7, height = 5)

sf_topo_stats %>%
  group_by(gamma) %>%
  summarise(sum_ss = sum(ss_total), .groups = 'drop') %>%
    ggplot(aes(x = gamma, y = sum_ss, fill = as.factor(gamma))) +
    geom_col() +
    geom_text(aes(label = sum_ss), vjust = -0.5, size = 7) +
    scale_y_continuous(labels = scales::label_number_si(accuracy = 0.1), limits = c(0,2500000)) +
    scale_x_continuous(breaks = c(2,2.5)) +
    scale_fill_brewer(palette = "Set1") +
    guides(fill = guide_legend(title = latex2exp::TeX("$\\gamma$"))) +
    labs(x = latex2exp::TeX("Scale-free exponent $\\gamma$"),
      title = "Total number of stable states", y = "") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "img/ss_dist_total_sum.png", dpi = "print", width = 7, height = 5)
