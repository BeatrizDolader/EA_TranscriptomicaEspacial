############################################################
#
# Identificación de HVGs en snRNAseq
#
############################################################

# 1. Paquetes necesarios y entorno Giotto
library(reticulate)
reticulate::use_condaenv("/clinicfs/userhomes/bdolader/.local/share/r-miniconda/envs/giotto_env", required = TRUE)
library(Giotto)

# 2. Cargar objeto Giotto
giotto_sn = Giotto::loadGiotto(path_to_folder = "/home/bdolader/TFM/01_snRNA/Giotto_sc")

# 3. Identificar HVGs
giotto_sn <- Giotto::calculateHVF(gobject = giotto_sn, 
                                  expression_values = "default_norm",
                                  method = "cov_loess",
                                  save_plot = TRUE)

# 4. Guardar objeto Giotto
saveGiotto(giotto_sn, dir = "/home/bdolader/TFM/01_snRNA", foldername = "Giotto_sc", overwrite = TRUE)
