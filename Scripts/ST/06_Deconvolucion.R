################################################################################
# 
# 06. Deconvolución
# 
# Este script realiza la deconvolución por el método DWLS con Giotto utilizando
# dos datasets: uno de snRNA y otro de ST; ambos cargados en un objeto Giotto.
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
# 2.1. Objetos Giotto (snRNA y ST)
giotto_ST = Giotto::loadGiotto(path_to_folder = "../../ST/Giotto")
giotto_sn = Giotto::loadGiotto(path_to_folder = "../../snRNA/Giotto_sn")

# 2.2. Genes marcadores
markers_file = "../../ST/gini_markers.rds"
cluster_markers = readRDS(markers_file)


# ==============================================================================
# 3. Deconvolución utilizando genes marcadores de grupo de Leiden
# ==============================================================================
# 3.1. Seleccionar los mejores marcadores
top_markers <- cluster_markers[, head(.SD, 50), by = "cluster"]$feats

# 3.2. Crear matriz de firmas celulares
DWLS_matrix <- makeSignMatrixDWLSfromMatrix(
                         matrix = getExpression(giotto_sn,
                         values = "default",
                         output = "matrix"),
                         cell_type = pDataDT(giotto_sn)$cell_type,
                         sign_gene = top_markers)

# 3.3. Método DWLS
giotto_ST <- Giotto::runDWLSDeconv(gobject = giotto_ST, 
                                   sign_matrix = DWLS_matrix,
                                   expression_values = "default",
                                   cluster_column = "leiden_0.5",
                                   n_cell = 10)


# ==============================================================================
# 4. Guardar objeto Giotto
# ==============================================================================
saveGiotto(giotto_ST, foldername = "Giotto", 
           dir = "../../ST/",
           overwrite = TRUE)


# ==============================================================================
# 5. Visualización de los resultados
# ==============================================================================
samples = unique(giotto_ST$list_ID)
metadata = pDataDT(giotto_ST)
colors = c(
  EX = "#377EB8",  
  INH = "#4DAF4A",  
  OPC = "#BC80BD", 
  PER = "#E6F5C9", 
  ASC = "#FFA500",  
  ODC = "#E7298A", 
  END = "#66C2A5", 
  MG = "#B22222",  
  SMC = "grey",  
  FBR = "#D55E00")   

spat_plots = list()
samples <- sort(samples)
for(sample in samples){
  cell_ids = metadata[list_ID == sample]$cell_ID
  giotto_obj = subsetGiotto(giotto_ST, cell_ids = cell_ids)
  spat = spatDeconvPlot(giotto_obj, 
                        show_image = FALSE,
                        radius = 90,
                        title = sample,
                        deconv_name = "DWLS",
                        cell_color_code = colors)
  spat_plots[[sample]] = spat}

combined_plot = wrap_plots(spat_plots[samples], ncol = 6)

ggsave(combined_plot, file = "../../ST/deconv.png", width = 16, height = 12)


# ==============================================================================
# 6. Deconvolución utilizando los genes marcadores de tipo celular
# ==============================================================================
cell_markers <- Giotto::findMarkers_one_vs_all(
  gobject = giotto_sn,
  method = "scran",
  expression_values = "raw",
  cluster_column = "cell_type",
  pval = 0.05,
  logFC = 0.5,
  min_feats = 20)

top <- cell_markers[, head(.SD, 50), by = "cluster"]
topgenes <- unique(top$feats)

DWLS_cell <- makeSignMatrixDWLSfromMatrix(
  matrix = getExpression(giotto_sn, values = "default", output = "matrix"),
  cell_type = pDataDT(giotto_sn)$cell_type,
  sign_gene = topgenes
)

st_genes = rownames(giotto_ST)
common_genes <- intersect(rownames(DWLS_cell), genes_st)
DWLS_cell_filtered <- DWLS_cell[genes_st[genes_st %in% rownames(DWLS_cell)], ]

giotto_ST <- Giotto::runDWLSDeconv(gobject = giotto_ST, 
                                   sign_matrix = DWLS_cell_filtered,
                                   name = "DWLS_cell",
                                   expression_values = "default",
                                   cluster_column = "leiden_0.5",
                                   n_cell = 10)
