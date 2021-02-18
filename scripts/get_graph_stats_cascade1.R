library(dplyr)
library(tibble)
library(readr)
library(igraph)

# CASCADE 1.0
edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_1.0.sif', delim = "\t", col_names = c('source', 'effect', 'target'), col_types = "ccc")

# data frame for igraph input
df = edge_tbl %>%
  rename(from = source, to = target) %>%
  relocate(effect, .after = to)

cascade_graph = igraph::graph_from_data_frame(df, directed = TRUE)

# get in-degrees, out-degrees and Clustering coefficients
in_degree = cascade_graph %>% igraph::degree(mode = 'in')
out_degree = cascade_graph %>% igraph::degree(mode = 'out')
all_degree = cascade_graph %>% igraph::degree(mode = 'all')
cluster_coef = cascade_graph %>% igraph::transitivity(type = 'local', isolates = 'zero')

# graph stats for CASCADE 1.0
cascade1_graph_stats = dplyr::bind_cols(node = V(cascade_graph)$name,
  in_degree = in_degree, out_degree = out_degree, all_degree = all_degree,
  cluster_coef = cluster_coef)

#' `cascade1_graph_stats` has for each node in the CASCADE 1.0 network
#' (a total of 77 nodes/rows) the number of inputs (`in_degree`), outputs
#' (`out_degree`), the sum of the two (`all_degree`) and the clustering
#' coefficient (`cluster_coef`)
saveRDS(cascade1_graph_stats, file = "data/cascade1_graph_stats.rds")
