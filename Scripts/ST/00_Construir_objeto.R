################################################################################
# 
# 00. Construcción del objeto Giotto
# 
# Este script crear un objeto Giotto a partir de datos de ST almacenados en un
# objeto Seurat.
#
################################################################################

# ==============================================================================
# 1. Paquetes 
# ==============================================================================
library(dplyr)
library(Giotto)
library(magick)
library(imager)


# ==============================================================================
# 2. Cargar datos descargados de GEO 
# ==============================================================================
seurat_object <- readRDS(paste0("../../Datos/",
                            "GSE233208_Human_visium_ADDS_seurat_processed.rds")


# ==============================================================================
# 3. Extraer información del objeto Seurat
# ==============================================================================
# 3.1. Matriz de conteos
raw_counts <- seurat_object[["Spatial"]]@counts
  
# 3.2. Metadatos
metadata = seurat_object@meta.data 
metadata$cell_ID <- colnames(raw_counts)
cell_metadata = metadata %>% dplyr::select(cell_ID, SampleID, Case.Year, 
                                           Case.Num, Diagnosis, Sex, Age, 
                                           Slide, seqbatch, Sample)

# 3.3. Coordenadas espaciales
coordenadas_list <- lapply(names(seurat_object@images), function(slice) {
  spatial_data <- seurat_object@images[[slice]]@coordinates
  data.frame(sample = slice,  
             cell_ID = rownames(spatial_data),
             sdimx = spatial_data$imagerow,
             sdimy = spatial_data$imagecol)})

coord_df <- do.call(rbind, coordenadas_list)

# 3.4. Imágenes histológicas
images <- list()

for (img_name in names(seurat_object@images)) {
  img_raster <- seurat_object@images[[img_name]]@image
  img_magick <- magick::image_read(img_raster)
  file_path <- file.path("../../ST/Images/", 
                         paste0(img_name, ".png"))

  # Guardar la imagen en formato .png
  magick::image_write(img_magick, file_path)
  
  # Leer la imagen
  img <- load.image(file_path)
  
  # Rotar la imagen 90 grados en sentido antihorario
  img_rotated <- imrotate(img, angle = -90)
  rotated_image_path <- file.path("../../ST/Images/", 
                                  paste0(img_name, "_rotated.png"))
  
  # Guardar la imagen rotada
  save.image(img_rotated, rotated_image_path)
  images[[img_name]] <- rotated_image_path}


# ==============================================================================
# 4. Filtrar los datos
# ==============================================================================
# 4. 1. Excluir las muestras de SD y réplicas biológicas
cell_metadata <- cell_metadata %>% mutate(Diagnosis = case_when(
                                          Diagnosis == "AD" ~ "EA",
                                          Diagnosis == "earlyAD" ~ "EA_temp",
                                          TRUE ~ Diagnosis))

groups = c("Control", "EA_temp", "EA")

to_discard = c("Dec_20_2021_Human1", "Dec_20_2021_Human2", "Dec_20_2021_Human3",
               "Dec_20_2021_Human4", "Dec_20_2021_Human5", "Dec_20_2021_Human6",
               "Dec_20_2021_Human7", "Dec_20_2021_Human8")

metadata_filtered <- cell_metadata %>%
  filter(Diagnosis %in% groups & !(Sample %in% to_discard))
  
counts_filtered <- raw_counts[, metadata_filtered$cell_ID]

coords_filtered <- coord_df %>%
  filter(cell_ID %in% metadata_filtered$cell_ID)

images_to_keep <- unique(coords_filtered$sample)

images_filtered <- images[names(images) %in% images_to_keep]


# ==============================================================================
# 5. Modificar el nombre de las muestras
# ==============================================================================
new_sampleid <- c("Oct_2021_1" = "Control_M_1",
                  "Oct_2021_2" = "EA_temp_H_1",
                  "Oct_2021_3" = "EA_M_1",
                  "Oct_2021_5" = "Control_H_1",
                  "Oct_2021_6" = "EA_temp_M_1",
                  "Oct_2021_7" = "EA_H_1",
                  "Nov_24_2021_VisiumHuman_1" = "Control_H_2",
                  "Nov_24_2021_VisiumHuman_2" = "EA_temp_M_2",
                  "Nov_24_2021_VisiumHuman_3" = "EA_M_2",
                  "Nov_24_2021_VisiumHuman_5" = "Control_M_2",
                  "Nov_24_2021_VisiumHuman_6" = "EA_temp_H_2",
                  "Nov_24_2021_VisiumHuman_7" = "EA_H_2",
                  "Nov_24_2021_VisiumHuman_9" = "Control_H_3",
                  "Nov_24_2021_VisiumHuman_10" = "EA_temp_M_3",
                  "Nov_24_2021_VisiumHuman_11" = "EA_H_3",
                  "Nov_24_2021_VisiumHuman_13" = "Control_M_3",
                  "Nov_24_2021_VisiumHuman_14" = "EA_temp_H_3",
                  "Nov_24_2021_VisiumHuman_15" = "EA_H_4",
                  "Dec_13_2021_Human1" = "Control_M_4",
                  "Dec_13_2021_Human2" = "EA_temp_H_4",
                  "Dec_13_2021_Human3" = "EA_M_3",
                  "Dec_13_2021_Human5" = "Control_H_4",
                  "Dec_13_2021_Human6" = "EA_temp_H_5",
                  "Dec_13_2021_Human7" = "EA_M_4")

metadata_filtered$SampleID <- new_sampleid[metadata_filtered$Sample]


# ==============================================================================
# 6. Construir objeto Giotto por muestra
# ==============================================================================
unique_samples <- unique(metadata_filtered$SampleID)
giotto_objects <- list()

for (sample in unique_samples) {
  sample_cells <- metadata_filtered[metadata_filtered$SampleID == sample, "cell_ID"]
  sample_counts <- counts_filtered[, sample_cells]
  sample_meta <- metadata_filtered[sample_cells, ]
  sample_coords <- coords_filtered[coords_filtered$cell_ID %in% sample_cells, ]
  sample_image <- images[names(images) %in% sample_coords$sample] 
  
  # Crear el objeto Giotto
  giotto_objects[[sample]] <- createGiottoObject(
    expression = sample_counts,
    spatial_locs = sample_coords[, c("cell_ID", "sdimx", "sdimy")],
    cell_metadata = sample_meta)
  
  # Crear objeto GiottoImage
  giotto_img = createGiottoImage(
    spatial_locs = sample_coords,
    mg_object = sample_image[[1]],
    name = names(sample_image))
  
  # Añadir imagen al objeto Giotto
  giotto_objects[[sample]] <- addGiottoImage(
    gobject = giotto_objects[[sample]],  
    images = list(giotto_img))
}


# ==============================================================================
# 7. Combinar las muestras en un único objeto Giotto
# ==============================================================================
giotto_combined <- Giotto::joinGiottoObjects(gobject_list = giotto_objects,
                                             gobject_names = names(giotto_objects),
                                             join_method = "shift",
                                             verbose = TRUE)

# ==============================================================================
# 8. Guardar objeto combinado
# ==============================================================================
saveGiotto(gobject = giotto_combined, 
           foldername = paste0("Giotto"),  
           dir = "../../ST/",
           method = "RDS",
           overwrite = TRUE,
           export_image = TRUE,
           image_filetype = "PNG",
           verbose = TRUE)

