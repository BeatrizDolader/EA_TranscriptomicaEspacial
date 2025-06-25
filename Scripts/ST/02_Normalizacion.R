################################################################################
# 
# 02. Normalización
# 
# Este script normaliza los conteos de un objeto Giotto mediante tamaño de 
# librería y el método de Pearson.
#
################################################################################


# ==============================================================================
# 1. Paquetes y entornos
# ==============================================================================
library(reticulate)
env_path = "/clinicfs/userhomes/bdolader/.local/share/r-miniconda/envs/giotto_env"
reticulate::use_condaenv(env_path, required = TRUE)
library(Giotto)
library(ggplot2)


# ==============================================================================
# 2. Cargar los datos
# ==============================================================================
giotto_ST = Giotto::loadGiotto(path_to_folder = "../../ST/Giotto")


# ==============================================================================
# 3. Normalización
# ==============================================================================
# 3.1. Por tamaño de librería
giotto_ST = Giotto::processExpression(giotto_ST, 
                                   param = normParam("default"),
                                   expression_values = "raw",
                                   name = "default")

# 3.2. Pearson
giotto_ST = Giotto::processExpression(giotto_ST, 
                                  param = normParam("pearson"),
                                  expression_values = "raw",
                                  name = "pearson")


# ==============================================================================
# 4. Guardar objeto Giotto
# ==============================================================================
saveGiotto(giotto_ST, foldername = "Giotto", 
           dir = "../../ST/",
           overwrite = TRUE)


# ==============================================================================
# 5. Visualización de los resultados
# ==============================================================================

colors <- c("Control" = "#1b9e77", "EA_temp" = "#d95f02", "EA" = "#7570b3")
fill_value_labels <- c("Control" = "Control", "EA_temp" = "EA temprana", 
                       "EA" = "EA tardía")

# 5.1. Conteos sin normalizar
metadata = pDataDT(giotto_ST)

metadata$Diagnosis <- factor(metadata$Diagnosis,
  levels = c("Control", "EA_temp", "EA"))

metadata$SampleID <- factor(metadata$SampleID,
  levels = unique(metadata$SampleID[order(metadata$Diagnosis)]))

box1 = ggplot(metadata, aes(x = SampleID, y = total_expr, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  labs(y = "Tamaño de librería", x = "Muestra", fill = "Grupo experimental") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 11)
  ) +
  scale_fill_manual(values = colors, labels = fill_value_labels)
ggsave(box1, file = "../../ST/raw_data.png", width = 16, height = 4, 
       dpi = 600)


# 5.2. Conteos normalizados
giotto_def = addCellStatistics(giotto_ST, expression_values = "default")
giotto_pearson = addCellStatistics(giotto_ST, expression_values = "pearson")

norm_list <- list(default = giotto_def, pearson = giotto_pearson)

for (norm_name in names(norm_list)) {
  giotto <- norm_list[[norm_name]]
  metadata_norm <- pDataDT(giotto)
  
  metadata_norm$Diagnosis <- factor(metadata_norm$Diagnosis,
                                    levels = c("Control", "EA_temp", "EA"))
  metadata_norm$SampleID <- factor(metadata_norm$SampleID,
      levels = unique(metadata_norm$SampleID[order(metadata_norm$Diagnosis)]))
  

  box <- ggplot(metadata_norm, aes(x = SampleID, y = total_expr, 
                                   fill = Diagnosis)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    labs(y = "Tamaño de librería", x = "Muestra", fill = "Grupo experimental") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.title = element_text(size = 11)
    ) +
    scale_fill_manual(values = colors, labels = fill_value_labels) 
  file <- paste0("../../ST/", norm_name, "_norm.png")
  ggsave(file, plot = box, width = 16, height = 4, dpi = 600)
}
