######################################
# Generate random topology SIF files #
######################################
library(BoolNet)
library(dplyr)
library(ggplot2)
library(readr)

# Define 2 useful functions first

# see BoolNet::generateRandomNKNetwork() for details about the input parameters
# This function returns an interaction tibble object (3 columns: source, effect
# and target), ready to be parsed into a .sif file
gen_random_nk_topology = function(n, k, topology, gamma = 2.5, approx_cutoff = 1000) {

  # define the vector of regulatory inputs (how many regulators per node in the network?)
  # based on the given input probability distribution [@Aldana2003]
  reg_inputs = switch(match.arg(arg = topology, choices = c("fixed", "homogeneous", "scale_free")),
    fixed = {
      if (length(k) == n) k
      else if (length(k) == 1) rep(k, n)
      else stop("k must have 1 or n element(s)!")
    },
    homogeneous = round(rpois(n, k)),
    scale_free = BoolNet:::rzeta(n, k, gamma = gamma, approx_cutoff = approx_cutoff),
    stop("`topology` must be one of `fixed`,`homogeneous` or `scale_free`")
  )

  # every node subject to regulation :)
  stopifnot(all(reg_inputs >= 1))

  # in case some value surpassed the total number of nodes
  reg_inputs[reg_inputs > n] = n

  # node names as `x` variables
  node_names = paste("x", seq_len(n), sep = "")
  names(reg_inputs) = node_names

  interaction_list = list()
  index = 1
  for (node in node_names) {
    regulator_num = reg_inputs[node]

    # get a random sample of regulators
    regulators = sample(x = node_names, size = regulator_num)

    for (regulator in regulators) {
      # source, effect, target
      df = dplyr::bind_cols(source = regulator,
        effect = sample(x = c("->", "-|"), size = 1), # `->` => activation, `-|` => inhibition
        target = node)
      interaction_list[[index]] = df
      index = index + 1
    }
  }

  stopifnot(length(interaction_list) == sum(reg_inputs))
  interaction_tbl = dplyr::bind_rows(interaction_list)

  return(interaction_tbl)
}

# write a .sif file with the results from `gen_random_nk_topology` function
write_sif = function(interactions_tbl, filename) {
  stopifnot(colnames(interactions_tbl) == c('source', 'effect', 'target'))
  readr::write_tsv(interactions_tbl, path = filename, col_names = FALSE)
}

##########################################
# Generate random homogeneous topologies #
##########################################
# n = 20 # #nodes
# k = 7 # mean number of regulators per node
#
# # for reproducibility
# set.seed(42)
# for (i in 1:100) {
#   res = gen_random_nk_topology(n = n, k = k, topology = "homogeneous")
#   write_sif(interactions_tbl = res, filename =
#       paste0("data/random_topology_files/homo_", n, "n_", k, "k_", i, ".sif"))
# }

#########################################
# Generate random scale free topologies #
#########################################
n = 50 # #nodes
k = 50 # max number of regulators per node

set.seed(42)
for (i in 1:100) {
  res = gen_random_nk_topology(n = n, k = k, topology = "scale_free") # gamma = 2.5 (default)
  write_sif(interactions_tbl = res, filename = paste0("data/random_topology_files/scale_free_", n, "n_", k, "k_gamma2_5_", i, ".sif"))
}

# gamma = 2 => hub nodes with even higher degree
set.seed(42)
for (i in 1:100) {
  res = gen_random_nk_topology(n = n, k = k, topology = "scale_free", gamma = 2)
  write_sif(interactions_tbl = res, filename = paste0("data/random_topology_files/scale_free_", n, "n_", k, "k_gamma2_", i, ".sif"))
}

# Other help functions
get_lo_stats = function(edge_tbl) {
  data_list = list()
  index = 1
  for (node in edge_tbl %>% distinct(target) %>% pull()) {
    effects = edge_tbl %>%
      filter(target == node) %>%
      pull(effect)
    if (length(unique(effects)) == 2) { # both activation and inhibition (node has link operator)
      num_reg = length(effects)
      num_act = sum(effects == "->")
      num_inh = sum(effects == "-|")
      data_list[[index]] = dplyr::bind_cols(lo_index = index, node = node, num_reg = num_reg,
        num_act = num_act, num_inh = num_inh)
      index = index + 1
    }
  }

  lo_stats = dplyr::bind_rows(data_list)
  #lo_stats %>% arrange(desc(num_reg))
}

print_stats = function() {
  for (topology_file in list.files(path = 'data/random_topology_files',
    full.names = TRUE, pattern = "scale")) {
    print(topology_file)
    edge_tbl = readr::read_delim(file = topology_file, delim = "\t",
      col_names = c('source', 'effect', 'target'), col_types = "ccc")
    lo_stats = get_lo_stats(edge_tbl)
    if (nrow(lo_stats) >= 20) {
      print(paste0("Number of link-operator equations: ", nrow(lo_stats) ,
        " with max #regulators: ", lo_stats %>% summarise(m = max(num_reg)) %>% pull()))
    }
  }
}

plot_distribution_fig = function(edge_tbl) {
  node_names = edge_tbl %>% pull(target) %>% unique()

  # calculate in-degrees
  in_degree = sapply(node_names, function(node) { edge_tbl %>% filter(target == node) %>% nrow() })

  # calculate out-degrees
  out_degree = sapply(node_names, function(node) { edge_tbl %>% filter(source == node) %>% nrow() })

  # degree distribution stats
  dd_stats = dplyr::bind_cols(node = node_names, in_degree = in_degree, out_degree = out_degree)

  dd_stats %>% group_by(in_degree) %>% tally() %>%
    ggplot(aes(x = in_degree, y = n)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    geom_smooth(aes(color = "red"), se = FALSE, show.legend = FALSE) +
    theme_classic() +
    labs(title = "In-Degree Distribution", x = "In Degree", y = "Number of Nodes")
}
