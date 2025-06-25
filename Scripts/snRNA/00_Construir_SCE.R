################################################################################
# 
# 00. Construir objeto SCE
# 
# Este script construye un objeto SingleCellExperiment con los datos de snRNA-seq
#
################################################################################

# ==============================================================================
# 1. Paquetes
# ==============================================================================
library(SingleCellExperiment)
library(Seurat)
library(dplyr)


# ==============================================================================
# 2. Leer datos
# ==============================================================================
seurat_object <- readRDS(paste0("../../Datos/",
                                "GSE233208_Human_snRNA-Seq_ADDS_integrated.rds")


# ==============================================================================
# 3. Extraer conteos y metadatos del objeto Seurat
# ==============================================================================
counts_matrix <- seurat_object@assays$RNA@counts
metadata = seurat_object@meta.data 
metadata = metadata %>% 
  select(cell_barcode, SampleID, Age, Sex, DX, NPDx1, Batch, cell_type, Tissue, 
         annotation, subtype)
metadata$barcode = rownames(metadata)


# ==============================================================================
# 4. Filtrar muestras
# ==============================================================================
group1 = "Alzheimer's disease"
group2 = "Control"

metadata_filtered = metadata[
  ((metadata$DX == group2) | (metadata$NPDx1 == group1)) & 
    (metadata$Tissue == "FCX"), ]

counts_filtered = counts_matrix[, colnames(counts_matrix) %in% 
                                  rownames(metadata_filtered)]


# ==============================================================================
# 5. Incluir nombre descriptivo a las muestras (Grupo_Sexo_Num)
# ==============================================================================
sample_metadata = metadata_filtered %>% distinct(SampleID, .keep_all = TRUE) %>% 
  select (SampleID, Age, Sex, DX, NPDx1)

sample_metadata <- sample_metadata %>%
  mutate(
    Grupo = ifelse(DX == "DSAD", "EA", "Control"),
    Sexo = ifelse(Sex == "F", "M", "H"),
    SampleTemp = paste0(Grupo, "_", Sexo)
  ) %>%
  group_by(SampleTemp) %>%
  mutate(
    Count = row_number(),
    Sample_name = paste0(SampleTemp, "_", Count)
  ) %>%
  ungroup() %>%
  select(-Grupo, -Sexo, -SampleTemp, -Count)

metadata_filtered <- metadata_filtered %>%
  left_join(sample_metadata %>% select(SampleID, Sample_name), by = "SampleID") %>%
  mutate(Grupo = ifelse(DX == "DSAD", "EA", "Control"))


# ==============================================================================
# 6. Construir objeto SingleCellExperiment
# ==============================================================================
genes = as.data.frame(rownames(counts_filtered))

sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_filtered), 
                                                 rowData = genes, 
                                                 colData = metadata_filtered)



# ==============================================================================
# 7. Guardar objeto SingleCellExperiment
# ==============================================================================
saveRDS(sce, file = "../../snRNA/sce.rds")