################################################################################
# 
# 04. Reducción de la dimensionalidad (PCA y UMAP)
# 
# Este script realiza la reducción de la dimensionalidad mediante PCA y su 
# posterior visualización con UMAP. Utiliza los tres conjuntos de HVGs
# identificados en el paso anterior.
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
giotto_ST = Giotto::loadGiotto(path_to_folder = "../../ST/Giotto")


# ==============================================================================
# 3. Reducción de dimensionalidad
# ==============================================================================
# 3.1. PCA
giotto_ST <- Giotto::runPCA(gobject = giotto_ST, 
                            expression_values = "pearson",
                            feats_to_use = "hvg_covgroup",
                            name = "pca_cov")

giotto_ST <- Giotto::runPCA(gobject = giotto_ST, 
                            expression_values = "pearson",
                            feats_to_use = "hvg_loess",
                            name = "pca_loess")

giotto_ST <- Giotto::runPCA(giotto_ST, 
                            expression_value = "pearson",
                            feats_to_use = "hvg_pearson",
                            name = "pca_pearson") 

# 3.2. Scree plot
Giotto::screePlot(giotto_ST, expression_values = "pearson", 
                  dim_reduction_name = "pca_cov", ncp = 50)

Giotto::screePlot(giotto_ST, expression_values = "pearson", 
                  dim_reduction_name = "pca_loess", ncp = 50)

Giotto::screePlot(giotto_ST, expression_values = "pearson", 
                  dim_reduction_name = "pca_pearson", ncp = 50)

# 3.3. UMAP
giotto_ST <- Giotto::runUMAP(giotto_ST, expression_values = "pearson", 
                             dim_reduction_name = "pca_cov", 
                             dimensions_to_use = 1:10, name = "umap_cov")

giotto_ST <- Giotto::runUMAP(giotto_ST, expression_values = "pearson", 
                             dim_reduction_name = "pca_loess", 
                             dimensions_to_use = 1:10, name = "umap_loess")

giotto_ST <- Giotto::runUMAP(giotto_ST, expression_values = "pearson", 
                             dim_reduction_name = "pca_pearson", 
                             dimensions_to_use = 1:10, name = "umap_pearson")


# ==============================================================================
# 4. Guardar objeto Giotto
# ==============================================================================
saveGiotto(giotto_ST, foldername = "Giotto", 
           dir = "../../ST/",
           overwrite = TRUE)


# ==============================================================================
# 5. Visualización de los resultados
# ==============================================================================
Giotto::plotUMAP(giotto_ST, 
                 dim_reduction_name = "umap_cov",
                 cell_color = "Diagnosis",
                 cell_color_code = c("#1b9e77", "#d95f02", "#7570b3"),
                 show_center_label = FALSE)

Giotto::plotUMAP(giotto_ST, 
                 dim_reduction_name = "umap_loess",
                 cell_color = "Diagnosis",
                 cell_color_code = c("#1b9e77", "#d95f02", "#7570b3"),
                 show_center_label = FALSE)

Giotto::plotUMAP(giotto_ST, 
                 dim_reduction_name = "umap_pearson",
                 cell_color = "Diagnosis",
                 cell_color_code = c("#1b9e77", "#d95f02", "#7570b3"),
                 show_center_label = FALSE)

dimPlot2D(giotto_ST, 
          dim_reduction_name = "umap_loess", label_size = 0,
          group_by = "Diagnosis", point_size = 0.7,
          cell_color = "Diagnosis",
          cell_color_code = c("Control" = "#1b9e77", "EA_temp" = "#d95f02", 
                              "EA" = "#7570b3"))