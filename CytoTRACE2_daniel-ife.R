#######
# CytoTRACE2 for daniel-ife
#######



set.seed(1234)
# Load data ---------------------------------------------------------------
daniel.ife <- readRDS("./rds/daniel_ife.rds")



# Import libraries --------------------------------------------------------
library("CytoTRACE2")



# Run CytoTRACE -----------------------------------------------------------
### Run main function ----
daniel.ife <- cytotrace2(
  input = daniel.ife, species = "mouse", is_seurat = TRUE, slot_type = "counts",
  full_model = TRUE, batch_size = 10000, smooth_batch_size = 1000,
  parallelize_models = TRUE, parallelize_smoothing = TRUE,
  ncores = NULL, max_pcs = 200, seed = 1234
)
daniel.ife@meta.data <- daniel.ife@meta.data %>% 
  dplyr::rename(
    cytotrace2.score = "CytoTRACE2_Score",  # absolute score from 0 (diff) to 1 (totipotent)
    cytotrace2.potency = "CytoTRACE2_Potency",  # potency assignment based on 6-binning of above
    cytotrace2.relative = "CytoTRACE2_Relative",  # relative score normalized to 0-1
    cytotrace2.preKNN_score = "preKNN_CytoTRACE2_Score",
    cytotrace2.preKNN_potency = "preKNN_CytoTRACE2_Potency"



#####
