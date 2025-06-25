################################################################################
# 
# 01. Control de calidad snRNA-seq
#  
# Este script realiza el control de calidad a los datos de snRNA-seq.
#
################################################################################

# ==============================================================================
# 1. Paquetes
# ==============================================================================
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scDblFinder)


# ==============================================================================
# 2. Cargar los datos
# ==============================================================================
sce = readRDS("../../snRNA/sce.rds")


# ==============================================================================
# 3. Control de calidad
# ==============================================================================
# 3.1. A nivel de célula
sce = scuttle::addPerCellQC(sce)
mt_genes <- grep("^MT-", rownames(sce), value = TRUE)
sce$ratio_mt = Matrix::colSums(counts(sce)[mt_genes, ]) / sce$total

# 3.2. A nivel de gen
sce = scuttle::addPerFeatureQC(sce) 

# 3.3. Identificación de dobletes
set.seed(123)
sce = scDblFinder::scDblFinder(sce, sample = "Sample_name", verbose = F)
table(sce$scDblFinder.class)

# 3.4. Filtrar genes que se expresan en muy pocas células
ncells_per_gene <- Matrix::rowSums(counts(sce) > 0)
x <- 10 
genes_to_keep <- ncells_per_gene >= x
sce_filtered <- sce[genes_to_keep, ]

# ==============================================================================
# 4. Guardar objeto SCE
# ==============================================================================
saveRDS(sce_filtered, file = "../../snRNA/sce.rds")


# ==============================================================================
# 5. Visualización del control de calidad
# ==============================================================================
  #  A nivel de núcleo
metadata = as.data.frame(colData(sce))

variables = c("total", "detected", "ratio_mt")
color_vars <- c("Grupo", "Batch")

fill_colors <- list(Grupo = c("Control" = "#1b9e77", "EA" = "#7570b3"), 
                    Batch = c("Batch1" = "#a6cee3", "Batch2" = "#b2df8a", 
                              "Batch3" = "#fb9a99", "Batch4" = "#fdbf6f", 
                              "Batch5" = "#cab2d6"))


y_labels <- c(total = "Tamaño de librería",
              detected = "Número de genes expresados",
              ratio_mt = "Ratio mitocondrial")

fill_labels <- c(Grupo = "Grupo experimental",
                 Batch = "Lote de secuenciación")

for (var in variables) {
  for (color_var in color_vars) {
    
    boxplot <- ggplot(metadata, aes(x = Sample_name, y = .data[[var]], 
                                    fill = .data[[color_var]])) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.8) + 
      labs(y = y_labels[[var]], x = "Muestra", fill = fill_labels[[color_var]]) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12)) +
      scale_fill_manual(values = fill_colors[[color_var]])
  
    file = paste0(var, "_", color_var, "_boxplot.png")
    ggsave(file, plot = boxplot, width = 10, height = 4, dpi = 300)}}


  # Dobletes
doublets_data <- metadata %>%
   group_by(Sample_name, scDblFinder.class) %>%
   tally(name = "count") %>%
   group_by(Sample_name) %>%
   mutate(prop = count / sum(count))

ggplot(doublets_data, aes(x = Sample_name, y = count, fill = scDblFinder.class)) +
  geom_bar(stat = "identity", position = "fill") +  
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(x = "Muestra", y = "Proporción (%)", fill = "Clase") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


