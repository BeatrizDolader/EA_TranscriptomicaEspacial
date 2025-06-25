################################################################################
# 
# 02. Normalización
# 
# Este script construye un objeto Giotto a partir de un objeto SCE y normaliza
# los conteos por tamaño de librería.
#
################################################################################


# ==============================================================================
# 1. Paquetes y entorno
# ==============================================================================
library(reticulate)
env_path = "/clinicfs/userhomes/bdolader/.local/share/r-miniconda/envs/giotto_env"
reticulate::use_condaenv(env_path, required = TRUE)
library(Giotto)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scDblFinder)

# ==============================================================================
# 2. Cargar los datos
# ==============================================================================
sce = readRDS("../../snRNA/sce.rds")


# ==============================================================================
# 3. Construir objeto Giotto
# ==============================================================================
counts = counts(sce)
metadata = colData(sce)

instructions = Giotto::createGiottoInstructions(
  save_dir = "/home/bdolader/TFM/01_snRNA",
  show_plot = TRUE,
  python_path = env_path)

giotto_sn <- Giotto::createGiottoObject(expression = counts,
                                        instructions = instructions)

giotto_sn <- addCellMetadata(giotto_sn,  new_metadata = metadata)
                            

# ==============================================================================
# 4. Normalización
# ==============================================================================
giotto_sn = Giotto::processExpression(giotto_sn, 
                                      param = normParam("default"),
                                      expression_values = "raw",
                                      name = "default")


# ==============================================================================
# 5. Guardar objeto Giotto 
# ==============================================================================
saveGiotto(giotto_sn, foldername = "Giotto_sn", 
           dir = "../../snRNA/",
           overwrite = TRUE)


# ==============================================================================
# 6. Visualización
# ==============================================================================
colors <- c("Control" = "#1b9e77", "EA" = "#7570b3")
metadata = pDataDT(giotto_sn)

raw = ggplot(metadata, aes(x = Sample_name, y = total, fill = Grupo)) +
             geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
             labs(y = "Tamaño de librería", x = "Muestra", 
                  fill = "Grupo experimental") +
             theme_minimal() +
             theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
              axis.text.y = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              legend.title = element_text(size = 11)
             ) +
             scale_fill_manual(values = colors)
ggsave(raw, file = "../../snRNA/raw_data.png", width = 10, height = 4, 
                    dpi = 600)

giotto_norm = addCellStatistics(giotto_sn, expression_values = "default")
metadata_norm <- pDataDT(giotto_norm)
  
normalized = ggplot(metadata_norm, aes(x = Sample_name, y = total_expr, fill = Grupo)) +
                    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
                    labs(y = "Tamaño de librería", x = "Muestra", 
                         fill = "Grupo experimental") +
                    theme_minimal() +
                    theme(
                      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                      axis.text.y = element_text(size = 14),
                      axis.title.x = element_text(size = 14),
                      axis.title.y = element_text(size = 14),
                      legend.title = element_text(size = 11)
                    ) +
                    scale_fill_manual(values = colors)
ggsave(normalized, file = "../../snRNA/normalized_data.png", width = 10, 
                         height = 4, dpi = 600)



