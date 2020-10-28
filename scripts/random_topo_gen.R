######################################
# Generate random topology SIF files #
######################################
library(BoolNet)
library(dplyr)
library(ggplot2)
library(readr)
source('scripts/help.R')

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
  readr::write_tsv(interactions_tbl, file = filename, col_names = FALSE)
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
for (i in 1:1000) {
  repeat {
    res = gen_random_nk_topology(n = n, k = k, topology = "scale_free") # gamma = 2.5 (default)
    lo_stats = get_lo_stats(res)
    if (nrow(lo_stats) < 19) {
      break # no more that 18 equations with link operators
    }
  }
  write_sif(interactions_tbl = res, filename = paste0("data/random_topology_files/scale_free_gamma2_5/scale_free_", n, "n_", k, "k_gamma2_5_", i, ".sif"))
}

# gamma = 2 => hub nodes with even higher degree
set.seed(42)
for (i in 1:100) {
  repeat {
    res = gen_random_nk_topology(n = n, k = k, topology = "scale_free", gamma = 2)
    lo_stats = get_lo_stats(res)
    if (nrow(lo_stats) < 19) {
      break # no more that 18 equations with link operators
    }
  }

  write_sif(interactions_tbl = res, filename = paste0("data/random_topology_files/scale_free_gamma2/scale_free_", n, "n_", k, "k_gamma2_", i, ".sif"))
}
