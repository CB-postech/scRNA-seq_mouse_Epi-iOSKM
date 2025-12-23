#######
# Script for running monocle3
#######



# Import libraries ----
library("monocle3")
library("SeuratWrappers")



# Load object ----
skin2.ife <- readRDS("./rds/skin2_ife.rds")



# Run monocle3 ----
### Make cds object ----
cds <- as.cell_data_set(skin2.ife)

### Modify dimension reductions ----
# the pipeline only supports the name "UMAP"
reducedDim(cds, "UMAP") <- reducedDim(cds, "TSNE.IFE")

### Cluster cells and learn graph ----
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

### Choose starting cell ----
skin2.ife$sc <- 
  skin2.ife$RNA$data["Krt14", ] *
  (
    max(skin2.ife$RNA$data["Krt10", ]) - skin2.ife$RNA$data["Krt10", ]
  )
sc_ <- which.max(skin2.ife$sc) %>% names()

### Order cells ----
cds <- order_cells(cds, root_cells = sc_)

### Plot result ----
plot_cells(
  cds = cds, 
  color_cells_by = "pseudotime",
  # color_cells_by = "cellType.epi",
  reduction = "UMAP",
  show_trajectory_graph = TRUE, trajectory_graph_segment_size = 1.5,
  trajectory_graph_color = "orangered1", 
  label_branch_points = FALSE, label_leaves = TRUE,
  label_cell_groups = FALSE,
  cell_size = 1, alpha = .25, cell_stroke = 1) + 
  NoLegend() + NoAxes() + ggtitle(NULL) +
  scale_color_viridis_c(option = "viridis") +
  theme(plot.margin = margin(.5, .5, .5, .5, "cm"))

### Save result to Seurat ----
# initial value
skin2.ife <- skin2.ife %>% 
  AddMetaData(
    object = .,
    metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
    col.name = "monocle3"
  )
# scale monocle3 values to a range 0-1
skin2.ife$monocle3.scaled <- (
  skin2.ife$monocle3 - min(skin2.ife$monocle3)
) / (
  max(skin2.ife$monocle3) - min(skin2.ife$monocle3)
)



#####
