library(obissdm)
library(parallel)
library(storr)
source("functions/model_utils.R")
source("functions/parallel_run.R")

# Load species list ----
all_sp_groups <- read.csv(list.files("data", pattern = "all_sp", full.names = T)[1])
# Include SDM groups based on the configuration file
all_sp_groups <- get_listbygroup(all_sp_groups)


# Create a storr for control ----
# This gives us more control over what was and what was not done
flux_ctrl_all <- storr_rds("sdmrun")

# Load habitat data ----
eco_info_all <- read.csv("data/species_ecoinfo.csv")



### Run in parallel
done_list <- storr::storr_rds("datastd")
done_list <- done_list$list()
#result <- mclapply(done_list, run_models, mc.cores = 5)
run_parallel(2, )


