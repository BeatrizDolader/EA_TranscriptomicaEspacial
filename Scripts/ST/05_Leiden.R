################################################################################
# 
# 05. Agrupamiento de Leiden
# 
# Este script aplica el algoritmo de Leiden para agrupar los spots en grupos 
# con perfiles de expresión génica similar.
#
################################################################################

# ==============================================================================
# 1. Paquetes y entornos
# ==============================================================================
library(reticulate)
env_path = "/clinicfs/userhomes/bdolader/.local/share/r-miniconda/envs/giotto_env"
reticulate::use_condaenv(env_path, required = TRUE)
library(Giotto)
library(patchwork)


# ==============================================================================
# 2. Cargar los datos
# ==============================================================================
giotto_ST = Giotto::loadGiotto(path_to_folder = "../../ST/Giotto")


# ==============================================================================
# 3. Red de vecinos más cercanos: sNN
# ==============================================================================
giotto_ST <- Giotto::createNearestNetwork(gobject = giotto_ST,
                                       dim_reduction_to_use = "pca",
                                       dim_reduction_name = "pca_loess",
                                       type = "sNN")


# ==============================================================================
# 4. Agrupamiento de Leiden
# ==============================================================================
res_values = c(0.3, 0.5, 0.7, 1)

for (value in res_values){
  name = paste0("leiden_", value)
  giotto_ST <- Giotto::doLeidenCluster(gobject = giotto_ST,
                                       name = name,
                                       nn_network_to_use = "sNN",
                                       network_name = "sNN.pca",
                                       resolution = value,
                                       n_iterations = -10)}


# ==============================================================================
# 5. Identificar marcadores génicos de grupo
# ==============================================================================
scran_markers <- Giotto::findMarkers_one_vs_all(gobject = giotto_ST,
                                                  method = "scran",
                                                  expression_values = "default",
                                                  cluster_column = "leiden_0.5",
                                                  pval = 0.01,
                                                  logFC = 0.5,
                                                  min_feats = 20)

gini_markers <- Giotto::findMarkers_one_vs_all(gobject = giotto_ST,
                                                method = "gini",
                                                expression_values = "default",
                                                cluster_column = "leiden_0.5")

# ==============================================================================
# 6. Guardar objeto Giotto y marcadores génicos
# ==============================================================================
saveGiotto(giotto_ST, foldername = "Giotto", 
           dir = "../../ST/",
           overwrite = TRUE)

saveRDS(scran_markers, file = "../../ST/scran_markers.rds")
saveRDS(gini_markers, file = "../../ST/gini_markers.rds")


# ==============================================================================
# 7. Visualización de los resultados
# ==============================================================================
# 7.1. Spat plot
palette_colors = c(
  "#FF0000", "#00BFFF", "#228833", "#FFD700", "#FF00FF", "#8A2BE2", "#FF6600", 
  "#B3DE69", "#FCCDE5", "#CC243C", "#44AAAA", "#FF6347", "#BBBBBB", "#1F77B4", 
  "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", 
  "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78")

names(palette_colors) = as.character(1:25)

samples = unique(giotto_ST$list_ID)
metadata = pDataDT(giotto_ST)
spat_plots = list()

for(sample in samples){
  cell_ids = metadata[list_ID == sample]$cell_ID
  giotto_obj = subsetGiotto(giotto_ST, cell_ids = cell_ids)
  spat = spatPlot(giotto_obj, 
                  cell_color = "leiden_0.5",
                  title = sample,
                  point_size = 1,
                  show_legend = FALSE,
                  cell_color_code = palette_colors)
  spat_plots[[sample]] = spat
  
}
combined_plot = wrap_plots(spat_plots[samples], ncol = 6)
ggsave(combined_plot, file = "spat_leiden_05.png", idth = 16, height = 12)

# 7.2. UMAP
umap_plot = Giotto::plotUMAP(giotto_ST, 
                             dim_reduction_name = "umap_loess",
                             cell_color = "leiden_0.5",
                             cell_color_code = palette_colors,
                             show_legend = TRUE)

ggsave(umap_plot, file = "umap_leiden_05.png", width = 16, height = 12)

# 7.3. Heatmap con los marcadores de grupo
topgenes_gini <- gini_markers[, head(.SD, 2), by = "cluster"]$feats
topgenes_scran <- scran_markers[, head(.SD, 2), by = "cluster"]$feats

plotMetaDataHeatmap(giotto_ST,
                    selected_feats = topgenes_gini,
                    expression_values = "default",
                    metadata_cols = "leiden_0.5",
                    custom_cluster_order = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 
                                             10, 11, 12, 13))
plotMetaDataHeatmap(giotto_ST,
                    selected_feats = topgenes_scran,
                    expression_values = "default",
                    metadata_cols = "leiden_0.5",
                    custom_cluster_order = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 
                                             10, 11, 12, 13))




