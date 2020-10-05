################################################################
# Calculate statistics (parameterization vs ss) for the
# link operator nodes in randomly-sampled CASCADE 2.0 models
#
# See Zenodo dataset: https://doi.org/10.5281/zenodo.3932382
################################################################

library(dplyr)
library(tibble)
library(emba)
library(usefun)

# Directory with all CASCADE 2.0 randomly-generated 1 stable-state models
data_dir = "/home/john/tmp/balance_paper/models_ss1"

ss_cascade2 = emba::get_stable_state_from_models_dir(models.dir = data_dir)
lo_cascade2 = emba::get_link_operators_from_models_dir(models.dir = data_dir)

# data check
stopifnot(rownames(ss_cascade2) == rownames(lo_cascade2))

# get link operator node statistics (number of regulators, activators, inhibitors)
edge_mat = emba::get_edges_from_topology_file(topology.file = "data/cascade_2_0.sif")
edge_tbl = edge_mat %>% as_tibble()

data_list = list()
index = 1
for (node in colnames(lo_cascade2)) {
  effects = edge_tbl %>%
    filter(target == node) %>%
    pull(regulation.effect)
  if (length(unique(effects)) == 2) { # both activation and inhibition (node has link operator)
    num_reg = length(effects)
    num_act = sum(effects == "activation")
    num_inh = sum(effects == "inhibition")
    data_list[[index]] = dplyr::bind_cols(lo_index = index, node = node, num_reg = num_reg,
      num_act = num_act, num_inh = num_inh)
    index = index + 1
  }
}

# `lo_stats` has info about the link operator nodes of CASCADE 2.0:
# Number of regulators (activators, inhibitors) and link operator equation index
# denoting the order we see the nodes when reading the logical rules
# (which is sequential since it was read from `lo_cascade2`)
lo_stats = dplyr::bind_rows(data_list)
lo_stats = lo_stats %>% mutate(across(where(is.double), as.integer))

# calculate `node_stats`
# It includes agreement measures between stable states and link operator assignments per node
data_list = list()
index = 1
for (node_name in lo_stats %>% pull(node)) {
  print(node_name)

  # get stable state data for the node (inhibited => 0, activated => 1)
  ss_data = ss_cascade2[, node_name]

  # get link operator data for the node (AND-NOT => 0, OR-NOT => 1)
  lo_data = lo_cascade2[, node_name]

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

  obs_prop_agreement = (and_not_0ss_agreement + or_not_1ss_agreement)/n
  p1 = ((or_not_1ss_agreement + and_not_1ss_disagreement)/n) * ((or_not_1ss_agreement + or_not_0ss_disagreement)/n)
  p0 = ((and_not_0ss_agreement + and_not_1ss_disagreement)/n) * ((and_not_0ss_agreement + or_not_0ss_disagreement)/n)
  random_prop_agreement = p1 + p0
  cohen_k = (obs_prop_agreement - random_prop_agreement)/(1 - random_prop_agreement)

  # calculate stats
  active_prop = sum(ss_data)/length(ss_data)
  inhibited_prop = 1 - active_prop
  or_not_prop = sum(lo_data)/length(lo_data)
  and_not_prop = 1 - or_not_prop

  data_list[[index]] = dplyr::bind_cols(
    obs_prop_agreement = obs_prop_agreement,
    cohen_k = cohen_k,
    and_not_0ss_agreement = and_not_0ss_agreement,
    or_not_1ss_agreement = or_not_1ss_agreement,
    and_not_1ss_disagreement = and_not_1ss_disagreement,
    or_not_0ss_disagreement = or_not_0ss_disagreement,
    active_prop = active_prop, inhibited_prop = inhibited_prop,
    or_not_prop = or_not_prop, and_not_prop = and_not_prop)
  index = index + 1
}

add_stats = dplyr::bind_rows(data_list)
node_stats = dplyr::bind_cols(lo_stats, add_stats)

# save result
# `node_stats` has a row for each node that has a link operator in its
# corresponding boolean equation (a total of 52 nodes/rows for CASCADE 2.0).
# The statistics included in each row have to do with the agreement between
# stable state activity and link operator parameterization for a random sample
# of CASCADE 2.0 models that have 1 stable state (a total of 20672 models)
saveRDS(node_stats, file = "data/node_stats_cascade2.rds")
