################################################################################
# 
# 01. Control de calidad
# 
# Este script realiza el control de calidad sobre un objeto Giotto con datos de
# ST.
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
library(patchwork)

# ==============================================================================
# 2. Cargar los datos
# ==============================================================================
giotto_ST = Giotto::loadGiotto(path_to_folder = "../../ST/Giotto")


# ==============================================================================
# 3. Control de calidad
# ==============================================================================
# 3. 1. A nivel de spot
giotto_ST = Giotto::addCellStatistics(giotto_ST, expression_values = "raw")

calculate_mt_ratio <- function(giotto_obj) {
  counts_matrix <- getExpression(giotto_obj, output = "matrix")
  mt_genes <- grep("^MT-", rownames(counts_matrix), value = TRUE)
  mt_counts <- counts_matrix[mt_genes, ]
  mt_counts_spot <- Matrix::colSums(mt_counts)
  metadata <- pDataDT(giotto_obj)
  ratio_mt <- mt_counts_spot / metadata$total_expr
  giotto_obj <- Giotto::addCellMetadata(
    gobject = giotto_obj,
    new_metadata = data.frame(ratio_mt = ratio_mt))
  
  return(giotto_obj)}

giotto_ST <- calculate_mt_ratio(giotto_ST)

# 3.3.1. Etiquetar los spots según su calidad
min_libsize <- 200
min_genes <- 100

metadata = pDataDT(giotto_ST)
total_expr_status <- ifelse(metadata$total_expr < min_libsize, "Descartado", 
"Aceptado")
nr_feats_status <- ifelse(metadata$nr_feats < min_genes, "Descartado", 
"Aceptado")
quality <- ifelse(total_expr_status == "Descartado" |
      nr_feats_status == "Descartado", "Baja calidad", "-")
  
new_pData <- data.frame(total_expr_status = as.factor(total_expr_status),
                        nr_feats_status = as.factor(nr_feats_status),
                        quality = as.factor(quality))
  
giotto_ST <- Giotto::addCellMetadata(gobject = giotto_ST, new_metadata = new_pData)
  

# 3.2. A nivel de gen
giotto_ST = addFeatStatistics(giotto_ST, expression_values = "raw")

# 3.2.1. Filtrar genes expresados en pocos spots
giotto_filtered = Giotto::filterGiotto(gobject = giotto_ST,
                                       expression_values = "raw",
                                       feat_det_in_min_cells = 3,
                                       min_det_feats_per_cell = 0)
# 3.2.2. Filtrar pseudogenes
total_genes = fDataDT(giotto_filtered)$feat_ID
pseudo <- "^(AC|AD|AF|AJ|AL|AP|BX|CR|CU|CY|FO|FP|KC|KF|L|U|Z)\\d+(\\.\\d+)?$"
genes_to_keep <- grep(pseudo, total_genes, value = TRUE, invert = TRUE)
giotto_filtered <- subsetGiotto(giotto_filtered, feat_ids = genes_to_keep)

  
# ==============================================================================
# 4. Guardar objet Giotto
# ==============================================================================
saveGiotto(giotto_filtered, foldername = "Giotto", 
           dir = "../../ST/", overwrite = TRUE)


# ==============================================================================
# 5. Visualización de los resultados
# ==============================================================================
# 5. 1. Violin y boxplot
metadata = pDataDT(giotto_ST)
variables <- c("total_expr", "nr_feats", "ratio_mt")
color_vars <- c("Diagnosis", "seqbatch")

fill_colors <- list(
  Diagnosis = c("Control" = "#1b9e77", "EA_temp" = "#d95f02", "EA" = "#7570b3"),
  seqbatch = c("Oct_2021" = "#a6cee3", "Nov_24_2021" = "#b2df8a", 
               "Dec_13_2021" = "#fb9a99"))

fill_value_labels <- list(
  Diagnosis = c("Control" = "Control", "EA_temp" = "EA temprana", 
                "EA" = "EA tardía"),
  seqbatch = c("Oct_2021" = "Batch 1", "Nov_24_2021" = "Batch 2", 
               "Dec_13_2021" = "Batch 3"))

y_labels <- c(total_expr = "Tamaño de librería",
              nr_feats = "Número de genes expresados",
              ratio_mt = "Ratio mitocondrial")

fill_labels <- c( Diagnosis = "Grupo experimental",
                  seqbatch = "Lote de secuenciación")

metadata$Diagnosis <- factor(metadata$Diagnosis,
                             levels = c("Control", "EA_temp", "EA"))

metadata$SampleID <- factor(metadata$SampleID,
  levels = unique(metadata$SampleID[order(metadata$Diagnosis)]))

for (var in variables) {
  for (fill_var in color_vars) {
  
    colors <- fill_colors[[fill_var]]
    fill_label <- fill_labels[[fill_var]]
    
    violin <- ggplot(metadata, aes_string(x = "SampleID", y = var, fill = fill_var)) +
      geom_violin(alpha = 0.7, scale = "width", trim = FALSE) +
      labs(y = y_labels[[var]], x = "Muestra", fill = fill_label) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 11)
      ) +
      scale_fill_manual(values = colors,
                        labels = fill_value_labels[[fill_var]])
    
    box <- ggplot(metadata, aes_string(x = "SampleID", y = var, fill = fill_var)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      labs(y = y_labels[[var]], x = "Muestra", fill = fill_label) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 11)
      ) +
      scale_fill_manual(values = colors,
                        labels = fill_value_labels[[fill_var]])
    
    plot <- violin / box
    file = paste0("../../ST/01_quality_control/", var, "_", fill_var, ".png")
    ggsave(file, plot = plot, width = 10, height = 8, dpi = 300)}}


# 5.2. Spots según su calidad
samples = unique(metadata$SampleID)
spat_plots = list()
for (sample in samples){
  cell_ids = metadata[list_ID == sample]$cell_ID
  giotto_obj = subsetGiotto(giotto_ST, cell_ids = cell_ids)
  
  p <- spatPlot(
    gobject = giotto_obj,
    show_image = FALSE,
    point_alpha = 1,
    cell_color = "quality",
    point_size = 1.5,
    cell_color_code = c("lightgreen", "red"),
    title = sample) +
    theme(legend.position = "none")
  
  spat_plots[[sample]] = p}

combined_plot = wrap_plots(spat_plots[samples], ncol = 6)
ggsave(combined_plot, file = "../../ST/01_quality_control/low_quality.png", 
       width = 16, height = 12)
