################################################################################
# 
# 07. Dominios espaciales
# 
# Este script crear una red espacial basada en la triangulación de Delaunay, 
# identifica los SVGs con BinSpect y determina los dominios espaciales bajo un
# modelo HMRF.
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
# 3. Red de Delaunay
# ==============================================================================
giotto_ST = createSpatialDelaunayNetwork(gobject = giotto_ST,
                                         method = "deldir")


# ==============================================================================
# 4. Genes espacialmente variables (SVGs)
# ==============================================================================
giotto_ST <- binSpect(gobject = giotto_ST, 
                      bin_method = "kmeans",
                      expression_values = "default",
                      return_gobject = TRUE)


# ==============================================================================
# 5. Dominios espaciales con HMRF
# ==============================================================================
init_k9 = initHMRF_V2(gobject = giotto_ST,
                      expression_values="default",
                      use_spatial_genes="binSpect",
                      k=9,
                      nstart=1000,
                      filter_method="none",
                      cl.method = "km",
                      resolution.cl = 1,
                      use_pca = FALSE, 
                      gene_list_from_top = 2500,
                      gene_samples = 500)

init_k10 = initHMRF_V2(gobject = giotto_ST,
                      expression_values="default",
                      use_spatial_genes="binSpect",
                      k=10,
                      nstart=1000,
                      filter_method="none",
                      cl.method = "km",
                      resolution.cl = 1,
                      use_pca = FALSE, 
                      gene_list_from_top = 2500,
                      gene_samples = 500)

res_k9 = doHMRF_V2(init_k9, betas=c(0, 2, 11))
res_k10 = doHMRF_V2(init_k10, betas = c(0, 2, 11))


# ==============================================================================
# 6. Guardar objeto Giotto
# ==============================================================================
saveGiotto(giotto_ST, foldername = "Giotto", 
           dir = "../../ST/",
           overwrite = TRUE)


# ==============================================================================
# 7. Visualización de los dominios
# ==============================================================================
giotto_ST = addHMRF_V2(giotto_ST, res_k9, name = "hmrf")
giotto_ST = addHMRF_V2(giotto_ST, res_k10, name = "hmrf")

samples = unique(giotto_ST$list_ID)
metadata = pDataDT(giotto_ST)
palette_colors = c(
  "#CC243C", "#00BFFF", "#228833", "#FFD700", "#FF00FF", "#8A2BE2", "#FF7F0E", 
  "#B3DE69", "#FCCDE5", "grey")
names(palette_colors) = as.character(1:10)

for(sample in samples){
  cell_ids = metadata[list_ID == sample]$cell_ID
  giotto_obj = subsetGiotto(giotto_ST, cell_ids = cell_ids)
  spat = spatPlot(giotto_obj, 
                  cell_color = "hmrf k=10 b=4.00",
                  title = sample,
                  point_size = 1,
                  show_legend = FALSE,
                  cell_color_code = palette_colors)
  spat_plots[[sample]] = spat
}

combined_plot = wrap_plots(spat_plots[samples], ncol = 6)

ggsave(combined_plot, file = "../../ST/k10_b4.png", width = 16, height = 12)
