#################################################
# Get link operator statistics data for all     #
# random scale free models generated via abmlog #
#################################################
library(dplyr)
library(readr)
library(stringr)
library(emba)
library(usefun)
source('scripts/help.R')

# directory with the `abmlog` results separated into 2 directories (per different gamma)
# and the scale free topology .sif files.
# See Zenodo dataset: https://doi.org/10.5281/zenodo.4392981
main_dir = "/media/disk/abmlog-scale-free/"

for (sf_dir in c("scale_free_gamma2_5", "scale_free_gamma2")) {
  gamma_param = ifelse(test = stringr::str_detect(string = sf_dir, pattern = "2_5"),
    yes  = 2.5, no = 2)
  print(paste0("Looking at results for gamma = ", gamma_param))

  res_dir = paste0(main_dir, sf_dir)
  topology_files = list.files(res_dir, pattern = ".sif", full.names = TRUE)

  # agreement stats (Link-operator vs Stable State)
  scale_free_models_stats_list = list()
  sf_index = 1

  # how many stable states had the models generated via abmlog from a single .sif file?
  abmlog_stats_list = list()
  ss_index = 1
  for (topology_file in topology_files) {
    # read the topology .sif file
    edge_tbl = readr::read_delim(file = topology_file, delim = "\t",
      col_names = c('source', 'effect', 'target'), col_types = "ccc")

    # get some link-operator stats
    lo_stats = get_lo_stats(edge_tbl)

    print(paste0("Topology file: ", basename(topology_file), " (", which(topology_files %in% topology_file), "/", length(topology_files), ") with ", nrow(lo_stats), " link operator equations"))

    # find abmlog result directory
    abmlog_dir_res = list.files(res_dir, pattern = paste0("results_",
      stringr::str_remove(string = basename(topology_file), pattern = ".sif"), "_"),
      full.names = TRUE)
    # check that we found one directory
    stopifnot(length(abmlog_dir_res) == 1)

    # get the model directories
    model_dirs = list.dirs(path = paste0(abmlog_dir_res, "/models"), recursive = FALSE)

    # get all stable state data from the abmlog-generated models
    print(paste0("Reading ", length(model_dirs), " model directories..."))
    ss_tbl_list = list()
    model_dir_index = 1
    for (model_dir in model_dirs) {
      #print(paste0("Reading directory: ", model_dir, " (", model_dir_index, ")"))
      ss_tbl_list[[model_dir_index]] = emba::get_stable_state_from_models_dir(models.dir = model_dir, all.ss = TRUE)
      model_dir_index = model_dir_index + 1
    }
    ss_tbl = dplyr::bind_rows(ss_tbl_list)

    # keep some stats
    abmlog_stats_list[[ss_index]] = dplyr::bind_cols(
      ss_total = nrow(ss_tbl),
      lo_nodes_num = nrow(lo_stats),
      max_num_reg = lo_stats %>% summarise(m = max(num_reg)) %>% pull(),
      gamma = gamma_param, file = basename(topology_file))
    ss_index = ss_index + 1

    # only continue if the parameterized abmlog-generated models had stable states
    if (nrow(ss_tbl) > 0) {
      print("Calculate agreement statistics (Link Operator vs Stable State Activity)")
      agreement_stats = calculate_agreement_stats(lo_stats, ss_tbl)
      node_stats = dplyr::bind_cols(lo_stats, agreement_stats)
      scale_free_models_stats_list[[basename(topology_file)]] = node_stats
      sf_index = sf_index + 1
    } else {
      print(paste0("No stable state data via abmlog for topology file: ", basename(topology_file)))
    }
  }

  abmlog_stats = dplyr::bind_rows(abmlog_stats_list)
  saveRDS(object = abmlog_stats, file = paste0("data/abmlog_stats_gamma", gamma_param, ".rds"))

  sf_stats = dplyr::bind_rows(scale_free_models_stats_list)
  saveRDS(object = sf_stats, file = paste0("data/sf_stats_gamma", gamma_param, ".rds"))
}
