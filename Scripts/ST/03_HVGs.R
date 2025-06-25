################################################################################
# 
# 03. Identificación de HVGs
# 
# Este script identifica genes altamente variables (HVGs) en un objeto Giotto
# mediante los tres métodos implementados por Giotto: covarianza por grupos, 
# regresión de LOESS y residuos de Pearson.
#
################################################################################


# ==============================================================================
# 1. Paquetes y entornos
# ==============================================================================
library(reticulate)
env_path = "/clinicfs/userhomes/bdolader/.local/share/r-miniconda/envs/giotto_env"
reticulate::use_condaenv(env_path, required = TRUE)
library(Giotto)
library(VennDiagram)


# ==============================================================================
# 2. Cargar los datos
# ==============================================================================
giotto_ST = Giotto::loadGiotto(path_to_folder = "../../ST/Giotto")


# ==============================================================================
# 3. Identificar HVGs
# ==============================================================================
# 3.1. Covarianza por grupos
giotto_ST <- Giotto::calculateHVF(gobject = giotto_ST, 
                                  expression_values = "default",
                                  method = "cov_group",
                                  HVFname = "hvg_covgroup", 
                                  save_plot = TRUE)

# 3.2. LOESS
giotto_ST <- Giotto::calculateHVF(gobject = giotto_ST, 
                                  expression_values = "default",
                                  method = "cov_loess",
                                  HVFname = "hvg_loess", 
                                  save_plot = TRUE)

# 3.3. Pearson
giotto_ST <- Giotto::calculateHVF(gobject = giotto_ST, 
                                  expression_values = "default",
                                  method = "var_p_resid",
                                  HVFname = paste0("hvg_pearson"), 
                                  save_plot = TRUE)


# ==============================================================================
# 4. Guardar objeto Giotto
# ==============================================================================
saveGiotto(giotto_ST, foldername = "Giotto", 
           dir = "../../ST/",
           overwrite = TRUE)


# ==============================================================================
# 5. Visualización de los resultados
# ==============================================================================
hvg_cols = c("hvg_covgroup", "hvg_loess", "hvg_pearson")
features_dt <- fDataDT(giotto_ST)
hvg_lists = list()

for (col in hvg_cols) {
  hvgs <- features_dt[features_dt[[col]] == "yes", ]$feat_ID
  hvg_lists[[col]] <- hvgs}

venn.diagram(x = hvg_lists, category.names = methods, filename = "venn_HVGs.png", 
             output = TRUE, imagetype = "png", height = 3000, width = 3000,
             resolution = 500, col = rep("black", length(methods)),  lwd = 2, 
             fill = c("#FCA192", "#AFDEA2", "#CFC982")[1:length(methods)],
             alpha = 0.5, cex = 1.8, cat.cex = 0, cat.pos = 0, cat.dist = 0.05,
             margin = 0.1)