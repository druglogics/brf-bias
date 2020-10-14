# some helpful functions
library(dplyr)
library(readr)
library(stringr)
library(usefun)
library(ggplot2)

get_topology_file_stats = function() {
  data_list = list()
  index = 1

  # read the scale free topology files with gamma = 2.5
  for (topology_file in list.files(path = 'data/random_topology_files/scale_free_gamma2_5',
    full.names = TRUE, pattern = "scale")) {
    # THIS IS HOW YOU READ A TAB-DELIMITED SIF FILE, RIGHT HERE!!!
    edge_tbl = readr::read_delim(file = topology_file, delim = "\t",
      col_names = c('source', 'effect', 'target'), col_types = "ccc")
    lo_stats = get_lo_stats(edge_tbl)

    data_list[[index]] = dplyr::bind_cols(lo_nodes_num = nrow(lo_stats), max_num_reg = lo_stats %>% summarise(m = max(num_reg)) %>% pull(), gamma = 2.5)
    index = index + 1
  }

  # read the scale free topology files with gamma = 2
  for (topology_file in list.files(path = 'data/random_topology_files/scale_free_gamma2',
    full.names = TRUE, pattern = "scale")) {
    # THIS IS HOW YOU READ A TAB-DELIMITED SIF FILE, RIGHT HERE!!!
    edge_tbl = readr::read_delim(file = topology_file, delim = "\t",
      col_names = c('source', 'effect', 'target'), col_types = "ccc")
    lo_stats = get_lo_stats(edge_tbl)

    data_list[[index]] = dplyr::bind_cols(lo_nodes_num = nrow(lo_stats), max_num_reg = lo_stats %>% summarise(m = max(num_reg)) %>% pull(), gamma = 2)
    index = index + 1
  }

  topo_file_stats = dplyr::bind_rows(data_list)
}

# `edge_tbl` is the result of reading a SIF file to a tibble (see above)
# The returned `lo_stats` has info about the link operator nodes of the model
# in the .sif file: Number of regulators (activators, inhibitors) and link
# operator equation index denoting the order we see the nodes when reading
# the logical rules
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
}

# calculate agreement measures between stable states and link operator assignments per node
# used in `get_scale_free_model_stats.R`
calculate_agreement_stats = function(lo_stats, ss_tbl) {
  # number of equations with link operator (`num_bits`)
  num_bits = nrow(lo_stats)

  # convert model numbers to unique link operator configuration in binary encoding
  # order of the models remains the same as in `ss_tbl`
  model_numbers = stringr::str_extract(ss_tbl %>% pull(model_name), pattern = "(\\d+$)") %>% as.integer()
  stopifnot(length(model_numbers) == nrow(ss_tbl))

  model_numbers_bin = sapply(model_numbers, usefun::dec_to_bin, num_bits)

  data_list = list()
  index = 1
  for (node_name in lo_stats %>% pull(node)) {
    #print(node_name)
    # get stable state data for the node
    ss_data = ss_tbl %>% pull(!!as.symbol(node_name))

    # which link operator equation has this node as target?
    # This index corresponds to the bit index in the binary representation
    node_index = lo_stats %>% filter(node == node_name) %>% pull(lo_index)

    # get link operator data for the node (AND-NOT => 0, OR-NOT => 1)
    lo_data = sapply(model_numbers_bin,
      function(num) {substr(num, node_index, node_index)}, USE.NAMES = FALSE) %>%
      as.numeric()

    # small data check
    stopifnot(length(ss_data) == length(lo_data))

    # calculate cohen's kappa (One rater assigns link-operator, another stable state!)
    add_ss = lo_data + ss_data
    sub_ss = lo_data - ss_data

    # AND-NOT (0) == Inhibited stable state (0)
    and_not_0ss_agreement = sum(add_ss == 0)
    # OR-NOT (1) == Active stable state (1)
    or_not_1ss_agreement = sum(add_ss == 2)
    # AND-NOT (0) != Active stable state (1)
    and_not_1ss_disagreement = sum(sub_ss == -1)
    # OR-NOT  (1) != Inhibited stable state (0)
    or_not_0ss_disagreement = sum(sub_ss == 1)

    n = and_not_0ss_agreement + or_not_1ss_agreement + and_not_1ss_disagreement + or_not_0ss_disagreement
    stopifnot(n == length(ss_data))

    percent_agreement = (and_not_0ss_agreement + or_not_1ss_agreement)/n
    p1 = ((or_not_1ss_agreement + and_not_1ss_disagreement)/n) * ((or_not_1ss_agreement + or_not_0ss_disagreement)/n)
    p0 = ((and_not_0ss_agreement + and_not_1ss_disagreement)/n) * ((and_not_0ss_agreement + or_not_0ss_disagreement)/n)
    random_prop_agreement = p1 + p0
    cohen_k = (percent_agreement - random_prop_agreement)/(1 - random_prop_agreement)

    # calculate stats
    active_prop = sum(ss_data)/length(ss_data)
    inhibited_prop = 1 - active_prop
    or_not_prop = sum(lo_data)/length(lo_data)
    and_not_prop = 1 - or_not_prop

    data_list[[index]] = dplyr::bind_cols(
      percent_agreement = percent_agreement,
      cohen_k = cohen_k,
      and_not_0ss_agreement = and_not_0ss_agreement,
      or_not_1ss_agreement = or_not_1ss_agreement,
      and_not_1ss_disagreement = and_not_1ss_disagreement,
      or_not_0ss_disagreement = or_not_0ss_disagreement,
      active_prop = active_prop, inhibited_prop = inhibited_prop,
      or_not_prop = or_not_prop, and_not_prop = and_not_prop)
    index = index + 1
  }

  agreement_stats = dplyr::bind_rows(data_list)
}

# in-degree distribution figure
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
