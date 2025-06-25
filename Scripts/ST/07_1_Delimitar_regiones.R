################################################################################
# 
# 07.1. Delimitar regiones espaciales
# 
# Este script crea una variable "Region" en el objeto Giotto que agrupa los 
# spots en sustancia gris o blanca según la comparación visual de los resultados
# de deconvolución y los dominios espaciales.
#
################################################################################

# ==============================================================================
# 1. Paquetes y entornos
# ==============================================================================
library(reticulate)
env_path = "/clinicfs/userhomes/bdolader/.local/share/r-miniconda/envs/giotto_env"
reticulate::use_condaenv(env_path, required = TRUE)
library(Giotto)


# ==============================================================================
# 2. Cargar los datos
# ==============================================================================
giotto_ST = Giotto::loadGiotto(path_to_folder = ".../../ST/Giotto")


# ==============================================================================
# 3. Unificar dominios en sustancia blanca y sustancia gris
# ==============================================================================
hmrf_clusters <- pDataDT(giotto_ST)$"hmrf k=10 b=4.00"

merged_clusters <- ifelse(hmrf_clusters %in% c("3", "4", "6", "10"), 
                          "Sustancia_blanca", 
                          "Sustancia_gris")

giotto_ST <- Giotto::addCellMetadata(giotto_ST,
                             new_metadata = merged_clusters,
                             vector_name = "Region")


samples = unique(giotto_ST$list_ID)
metadata = pDataDT(giotto_ST)

for(sample in samples){
  cell_ids = metadata[list_ID == sample]$cell_ID
  giotto_obj = subsetGiotto(giotto_ST, cell_ids = cell_ids)
   spat = spatPlot(giotto_obj, 
                  cell_color = "Region", 
                  title = sample,
                  point_size = 2.1,
                  cell_color_code = c("white", "grey"))
  spat_plots[[sample]] = spat
}
combined_plot = wrap_plots(spat_plots[samples], ncol = 6)
ggsave(combined_plot, file = "../../ST/matter_comp.png", idth = 16, height = 12)


# ==============================================================================
# 4. Guardar objeto Giotto 
# ==============================================================================
saveGiotto(giotto_ST, foldername = "Giotto", 
           dir = "../../ST/",
           overwrite = TRUE)