#################################################
# Get 1 stable state model data for CASCADE 1.0 #
#################################################
library(dplyr)
library(emba)

# point to the `models` directory in the Zenodo dataset (https://doi.org/10.5281/zenodo.4022783)
data_dir = "/media/disk/abmlog/abmlog_cascade_1.0_models_fixpoints/models"
model_dirs = list.dirs(path = data_dir, recursive = FALSE)

ss_data = list()
index = 1
for (model_dir in model_dirs) {
  print(paste0("Reading directory: ", model_dir, " (", index, ")"))
  ss_data[[index]] = emba::get_stable_state_from_models_dir(models.dir = model_dir)
  index = index + 1
}

res = dplyr::bind_rows(ss_data)

# `res` (`ss_data`) is a data.frame with rows the boolean models that had one
# stable state (~3 million) and columns the nodes of the CASCADE 1.0 network (77).
# Each value in the data.frame is either 0 (inhibited state) or 1 (active state).
saveRDS(res, file = "data/ss_data.rds")
