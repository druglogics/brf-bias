library(dplyr)
library(tibble)
library(readr)

# CASCADE 1.0
edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_1.0.sif', delim = "\t", col_names = c('source', 'effect', 'target'), col_types = "ccc")

sources = edge_tbl %>% pull(source)
targets = edge_tbl %>% pull(target)
nodes = unique(c(sources, targets))

# calculate in-degrees
in_degree = sapply(nodes, function(node) {
  edge_tbl %>% filter(target == node) %>% nrow()
})

# calculate out-degrees
out_degree = sapply(nodes, function(node) {
  edge_tbl %>% filter(source == node) %>% nrow()
})

# degree distribution stats for CASCADE 1.0
dd_stats = dplyr::bind_cols(node = nodes, in_degree = in_degree, out_degree = out_degree)

# `dd_stats` has for each node in the CASCADE 1.0 network (a total of 77 nodes/rows)
# the number of inputs (`in_degree`) and outputs (`out_degree`)
saveRDS(dd_stats, file = "data/dd_stats.rds")
