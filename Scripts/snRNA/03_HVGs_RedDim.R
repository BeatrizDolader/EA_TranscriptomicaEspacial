################################################################################
# 
# 03. Identificación de HVGs y reducción de dimensionalidad
# 
# Este script identifica genes altamente variables y reduce la dimensionalidad
# de los datos mediante PCA y UMAP.
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
giotto_sn = Giotto::loadGiotto(path_to_folder = "../../snRNA/Giotto_sn")


# ==============================================================================
# 3. Identificar HVGs
# ==============================================================================
giotto_sn <- Giotto::calculateHVF(gobject = giotto_sn, 
                                  expression_values = "default",
                                  method = "cov_loess",
                                  save_plot = TRUE)


# ==============================================================================
# 4. PCA y UMAP
# ==============================================================================
giotto_sn <- Giotto::runPCA(gobject = giotto_sn, 
                            expression_values = "default",
                            feats_to_use = "hvf")

Giotto::screePlot(giotto_sn, expression_values = "default", 
                  dim_reduction_name = "pca", ncp = 50)

giotto_sn <- Giotto::runUMAP(giotto_sn, expression_values = "default", 
                             dim_reduction_name = "pca_loess", 
                             dimensions_to_use = 1:10)

Giotto::plotUMAP(giotto_sn, 
                 dim_reduction_name = "umap",
                 cell_color = "cell_type")


# ==============================================================================
# 5. Guardar objeto Giotto
# ==============================================================================
saveGiotto(giotto_sn, dir = "../../snRNA", 
           foldername = "Giotto_sn", overwrite = TRUE)
