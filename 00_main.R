#############################################################################
# THIS IS THE MAIN SCRIPT, RUNNING ALL OTHER SUB_SCRIPTS
# NO NEED TO RUN THE SCRIPTS INDIVIDUALLY
############################################################################
## set-up


#empty environment to make sure no old variables are msitakenly used
rm(list=ls())
myseed<-100
set.seed(myseed)

# sets the working directory to the current location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## for repeat testing:
# # move old output to archive
# archive_folder <- paste0("archive/output_", format(Sys.time(), '%Y-%m-%d_%H-%M-%S'))
# dir.create(archive_folder)
# file.copy(from="output",
#           to=archive_folder, 
#           overwrite = TRUE, recursive = TRUE, 
#           copy.mode = TRUE)
# 
# ## empty output folder
# # remove old output folder
# unlink("output", recursive = TRUE)
# # create new output folder
dir.create("output")


# choice of scenario:
scenario <- "scen_1"

scen_1<-"Complete model with all household interventions and all treatment options"
scen_2<-"Scenario: all household interventions, no treatment"
scen_3<-"Scenario: all household interventions, radical cure treatment for everyone"
scen_4<-"Scenario: all household interventions, standard treatment for everyone"

# run all scripts in the right sequence
source("01_plot_MAP.R")
source("02_set_up_sim_code.R")
source(paste0('03_model_params_', scenario, '.R'))
source("04_multi-STABS_model.R")
source("05_plots.R")

# save workspace data if needed:
# save.image(paste0("output/", scenario, ".RData"))

