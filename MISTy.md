## MISTy



```bash
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02663-5

https://www.nature.com/articles/s41586-022-05060-x

https://github.com/saezlab/visium_heart/tree/master/st_snRNAseq/05_colocalization

https://github.com/saezlab/visium_heart/blob/master/st_snRNAseq/utils/misty_utilities.R

https://www.nature.com/articles/s41588-024-01914-4

https://saezlab.github.io/mistyR/articles/mistyR.html

https://saezlab.github.io/mistyR/articles/FunctionalPipelinePathwayActivityLigands.html

https://saezlab.github.io/mistyR/articles/FunctionalPipelinePathwaySpecific.html

https://saezlab.github.io/mistyR/articles/MistyRStructuralAnalysisPipelineC2L.html


```



```bash
# Set Bioconductor mirror and install necessary package
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install("mistyR")

# Load required libraries
library(mistyR)        # MISTy package
library(future)        # Parallel processing
library(dplyr)         # Data manipulation
library(purrr)         # Functional programming
library(distances)     # Distance calculations
library(ggplot2)       # Plotting
setwd("")  # Set working directory

# Set up parallel processing
plan(multisession)

# Load synthetic data and visualize
data("synthetic")
View(synthetic[["synthetic1"]])

ggplot(synthetic[[1]], aes(x = col, y = row, color = type)) +
  geom_point(shape = 15, size = 0.7) +
  scale_color_manual(values = c("#e9eed3", "#dcc38d", "#c9e2ad", "#a6bab6")) +
  theme_void()

# Create MISTy views and perform summary
expr <- synthetic[[1]] %>% select(-c(row, col, type, starts_with("lig")))
misty.intra <- create_initial_view(expr)
summary(misty.intra)
summary(misty.intra$intraview)

# Define positions and add paraviews
pos <- synthetic[[1]] %>% select(row, col)
misty.views <- misty.intra %>% add_paraview(pos, l = 10)
summary(misty.views)

# Create nearest neighbor view and add to MISTy views
neighbors <- nearest_neighbor_search(distances(as.matrix(pos)), k = 11)[-1, ]
nnexpr <- seq_len(nrow(expr)) %>%
  map_dfr(~ expr %>% slice(neighbors[, .x]) %>% colMeans())
nn.view <- create_view("nearest", nnexpr, "nn")
extended.views <- misty.views %>% add_views(nn.view)
summary(extended.views)

# Remove views and rerun MISTy
extended.views %>% remove_views("nearest") %>% summary()
extended.views %>% remove_views("intraview") %>% summary()
misty.views %>% run_misty()

# Run MISTy pipeline for synthetic samples
result.folders <- synthetic %>% imap_chr(function(sample, name) {
  sample.expr <- sample %>% select(-c(row, col, type, starts_with("lig")))
  sample.pos <- sample %>% select(row, col)
  create_initial_view(sample.expr) %>%
    add_paraview(sample.pos, l = 10) %>%
    run_misty(results.folder = paste0("results", .Platform$file.sep, name))
})
result.folders

# Process and summarize results
misty.results <- collect_results(result.folders)
summary(misty.results)
str(misty.results)

# Plotting result statistics and contributions
misty.results %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")
  
gain.R2 = multi.R2 - intra.R2
gain.RMSE = 100 * (intra.RMSE - multi.RMSE) / intra.RMSE

misty.results$improvements %>%
  filter(measure == "p.R2") %>%
  group_by(target) %>%
  summarize(mean.p = mean(value)) %>%
  arrange(mean.p)

misty.results %>% plot_view_contributions()
misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.8)
misty.results %>% plot_interaction_heatmap(view = "para.10", cutoff = 0.5)
misty.results %>% plot_contrast_heatmap("intra", "para.10", cutoff = 0.5)
misty.results %>% plot_interaction_communities("intra")
misty.results %>% plot_interaction_communities("para.10", cutoff = 0.5)

```



```bash
# MISTy
library(mistyR)
library(future)
#Seurat
library(Seurat)
library(SeuratObject)
# Data manipulation
library(tidyverse)
# Distances
library(distances)
setwd("10X_Visium_ACH005/ACH005/")
seurat_vs <- readRDS("ACH005.rds")


# Extract cell composition and location data
composition <- as_tibble(t(seurat_vs[["c2l_props"]]$data))
geometry <- GetTissueCoordinates(seurat_vs, cols = c("imagerow", "imagecol"), scale = NULL)

# Spatial plot for visual inspection
SpatialPlot(seurat_vs, keep.scale = NULL, alpha = 0)

# Cell type proportions and spatial feature plot
DefaultAssay(seurat_vs) <- "c2l_props"
SpatialFeaturePlot(seurat_vs, keep.scale = NULL, features = "CM")

# Calculate neighborhood radius
geom_dist <- as.matrix(distances(geometry))
dist_nn <- apply(geom_dist, 1, function(x) sort(x)[2])
paraview_radius <- ceiling(mean(dist_nn + sd(dist_nn)))

# Create initial view and add paraview with Gaussian kernel
heart_views <- create_initial_view(composition) %>%
  add_paraview(geometry, l = paraview_radius, family = "gaussian")

# Run MISTy and collect results
run_misty(heart_views, "result/vignette_structural_pipeline")
misty_results <- collect_results("result/vignette_structural_pipeline")

# Downstream analysis and plotting
misty_results %>%
  plot_improvement_stats("multi.R2") %>%
  plot_improvement_stats("gain.R2")

misty_results %>% plot_interaction_heatmap(view = "intra", clean = TRUE)

# Filter for specific predictors and plot
misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "CM") %>%
  arrange(-Importance)

SpatialFeaturePlot(seurat_vs, keep.scale = NULL, features = c("Fib", "CM"), image.alpha = 0)

# Additional interaction and community plots
misty_results %>% plot_interaction_heatmap(
  view = "para.126",
  clean = TRUE,
  trim = 1.75,
  trim.measure = "gain.R2",
  cutoff = 0.5
)

SpatialFeaturePlot(seurat_vs, keep.scale = NULL, features = c("prolif", "Adipo"), image.alpha = 0)
misty_results %>% plot_view_contributions()
misty_results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.8)
misty_results %>% plot_interaction_heatmap(view = "para.126", cutoff = 0.5)
misty_results %>% plot_contrast_heatmap("intra", "para.126", cutoff = 0.5)
misty_results %>% plot_interaction_communities("intra")
misty_results %>% plot_interaction_communities("para.126", cutoff = 0.5)

```





