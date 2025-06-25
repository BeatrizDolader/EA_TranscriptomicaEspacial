################################################################################
# 
# 08. Análisis de expresión diferencial con Scran.
# 
# Este script analiza las diferencias en la expresión en la sustancia gris en 3 
# comparaciones:
#       1) Mujeres control vs Hombres control
#       2) Mujeres EA temprana vs. Hombres EA temprana
#       3) Mujeres EA tardía vs. Hombres EA tardía
#
################################################################################


# ==============================================================================
# 1. Paquetes y entornos
# ==============================================================================
library(reticulate)
library(Giotto)
env_path = "/clinicfs/userhomes/bdolader/.local/share/r-miniconda/envs/giotto_env"
reticulate::use_condaenv(env_path, required = TRUE)
library(ggplot2)
library(VennDiagram)


# ==============================================================================
# 2. Cargar los datos
# ==============================================================================
giotto_ST = Giotto::loadGiotto(path_to_folder = "../../Giotto")


# ==============================================================================
# 3. Seleccionar los datos de sustancia gris
# ==============================================================================
region = "Sustancia_gris"
cell_ids = metadata[Region == region]$cell_ID
giotto_GM = subsetGiotto(giotto_ST, cell_ids = cell_ids)
metadata_GM = pDataDT(giotto_GM)


# ==============================================================================
# 4. Análisis de expresión diferencial
# ==============================================================================
# 4.1. Comparación 1: Mujeres control vs. Hombres control
cell_control = metadata_GM[Diagnosis == "Control"]$cell_ID
giotto_control = subsetGiotto(giotto_GM, cell_ids = cell_control)

markers_control = Giotto::findMarkers(giotto_control,
                                      expression_values = "default",
                                      cluster_column = "Sex",
                                      method = "scran")

sign_control = markers_control[[1]][FDR < 0.05 & abs(summary.logFC) > 0.5]
up_in_females <- sign_control[summary.logFC > 0]
up_in_males <- sign_control[summary.logFC < 0]

# 4.2. Comparación 2: Mujeres EA_temp vs. Hombres EA_temp
cell_EAtemp = metadata_GM[Diagnosis == "EA_temp"]$cell_ID
giotto_EAtemp =  subsetGiotto(giotto_GM, cell_ids = cell_EAtemp)

markers_EAtemp = Giotto::findMarkers(giotto_EAtemp,
                                      expression_values = "default",
                                      cluster_column = "Sex",
                                      method = "scran")

sign_EAtemp = markers_EAtemp[[1]][FDR < 0.05 & abs(summary.logFC) > 0.5]
up_EAtemp_F <- sign_EAtemp[summary.logFC > 0]
up_EAtemp_M <- sign_EAtemp[summary.logFC < 0]

# 4.3. Comparación 3: Mujeres EA tardía vs. Hombres EA tardía
cell_EA = metadata_GM[Diagnosis == "EA"]$cell_ID
giotto_EA =  subsetGiotto(giotto_GM, cell_ids = cell_EA)

markers_EA = Giotto::findMarkers(giotto_EA,
                                     expression_values = "default",
                                     cluster_column = "Sex",
                                     method = "scran")

sign_EA = markers_EA[[1]][FDR < 0.05 & abs(summary.logFC) > 0.5]
up_EA_F <- sign_EA[summary.logFC > 0]
up_EA_M <- sign_EA[summary.logFC < 0]


# ==============================================================================
# 5. Visualización
# ==============================================================================
# 5.1. Venn Diagram
# Genes sobreexpresados en mujeres
png("../../ST/venn_up_females.png", width = 1600, height = 1600, res = 300)
venn_f <- venn.diagram(
  x = list(
    Control = up_in_females$feats,
    EA_temp = up_EAtemp_F$feats,
    EA = up_EA_F$feats),
  filename = NULL,
  output = TRUE,
  fill = c("#A8E6A1", "#FFD59E", "#D2B3E5"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  category.names = c("","",""),
  main = "Genes sobreexpresados en mujeres")
grid::grid.draw(venn_f)
dev.off()

# Genes sobreexpresados en hombres
png("../../ST/venn_up_males.png", width = 1600, height = 1600, res = 300)
venn_m <- venn.diagram(
  x = list(
    Control = up_in_males$feats,
    EA_temp = up_EAtemp_M$feats,
    EA = up_EA_M$feats),
  filename = NULL,
  output = TRUE,
  fill = c("#A8E6A1", "#FFD59E", "#D2B3E5"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  category.names = c("","",""),
  main = "Genes sobreexpresados en hombres")
grid::grid.draw(venn_m)
dev.off()

# 5.2. Volcano plot
comparisons <- list(
  Control = markers_control[[1]],
  EA_temp = markers_EAtemp[[1]],
  EA = markers_EA[[1]])

color_map <- c(
  "Up in females" = "red",
  "Up in males" = "blue",
  "Not significant" = "gray80")

for (name in names(comparisons)) {
  volcano_data <- comparisons[[name]]

  volcano_data$significance <- with(volcano_data, ifelse(
    FDR < 0.05 & summary.logFC > 0.5, "Up in females",
    ifelse(FDR < 0.05 & summary.logFC < -0.5, "Up in males", "No significativos")))
  
  volcano_data$FDR[volcano_data$FDR == 0] <- 1e-320
  top_up_females <- volcano_data[volcano_data$significance == "Up in females", ]
  top_up_females <- top_up_females[order(-top_up_females$summary.logFC), ][1:10, ]
  top_up_males <- volcano_data[volcano_data$significance == "Up in males", ]
  top_up_males <- top_up_males[order(top_up_males$summary.logFC), ][1:10, ]
  top_genes <- rbind(top_up_females, top_up_males)
  
  p <- ggplot(volcano_data, aes(x = summary.logFC, y = -log10(FDR))) +
    geom_point(aes(color = significance), alpha = 0.7) +
    scale_color_manual(values = color_map) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_text_repel(
      data = top_genes,
      aes(label = feats),
      size = 5,         
      max.overlaps = 20
    ) +
    labs(
      x = "log2 Fold Change",
      y = "-log10(FDR)",
      color = "Grupo"
    ) +
    coord_cartesian(ylim = c(0, 350)) +
    theme_minimal() +
    theme(legend.position = "none") 
    
  ggsave(
    filename = paste0("volcano_", name, ".png"),
    plot = p,
    width =8, height = 8, dpi = 300)
}