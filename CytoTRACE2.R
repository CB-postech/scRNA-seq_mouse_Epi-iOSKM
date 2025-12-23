#######
# Script used in CytoTRACE analysis
#######


set.seed(1234)
# Import libraries --------------------------------------------------------
library("CytoTRACE2")



# Load data ---------------------------------------------------------------
skin1.ife <- readRDS("./rds/2024/0415/skin1_ife.rds")
skin2.ife <- readRDS("./rds/2024/0415/skin2_ife.rds")



# Run CytoTRACE -----------------------------------------------------------
# skin1.ife
skin1.ife <- cytotrace2(
  input = skin1.ife, species = "mouse", is_seurat = TRUE, slot_type = "counts",
  full_model = TRUE, batch_size = 10000, smooth_batch_size = 1000,
  parallelize_models = TRUE, parallelize_smoothing = TRUE,
  ncores = NULL, max_pcs = 200, seed = 1234
  )
skin1.ife@meta.data <- skin1.ife@meta.data %>% 
  dplyr::rename(
    cytotrace2.score = "CytoTRACE2_Score",  # absolute score from 0 (diff) to 1 (totipotent)
    cytotrace2.potency = "CytoTRACE2_Potency",  # potency assignment based on 6-binning of above
    cytotrace2.relative = "CytoTRACE2_Relative",  # relative score normalized to 0-1
    cytotrace2.preKNN_score = "preKNN_CytoTRACE2_Score",
    cytotrace2.preKNN_potency = "preKNN_CytoTRACE2_Potency"
  )
# skin2.ife
skin2.ife <- cytotrace2(
  input = skin2.ife, species = "mouse", is_seurat = TRUE, slot_type = "counts",
  full_model = TRUE, batch_size = 10000, smooth_batch_size = 1000,
  parallelize_models = TRUE, parallelize_smoothing = TRUE,
  ncores = NULL, max_pcs = 200, seed = 1234
)
skin2.ife@meta.data <- skin2.ife@meta.data %>% 
  dplyr::rename(
    cytotrace2.score = "CytoTRACE2_Score",  # absolute score from 0 (diff) to 1 (totipotent)
    cytotrace2.potency = "CytoTRACE2_Potency",  # potency assignment based on 6-binning of above
    cytotrace2.relative = "CytoTRACE2_Relative",  # relative score normalized to 0-1
    cytotrace2.preKNN_score = "preKNN_CytoTRACE2_Score",
    cytotrace2.preKNN_potency = "preKNN_CytoTRACE2_Potency"
  )



#####
