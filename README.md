# Transcriptómica con resolución espacial aplicada al estudio de las diferencias de sexo en la enfermedad de Alzheimer

## Trabajo Fin de Máster en Bioinformática - Universidad de Valencia

Este repositorio contiene el código desarrollado para el análisis de un conjunto de datos combinados de transcriptómica espacial (ST) y snRNA-seq, publicados por [Miyoshi et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11631771/), con el objetivo de estudiar las diferencias relacionadas con el sexo en la enfermedad de Alzheimer.

El análisis se organiza en dos secciones principales, cada una con su propio directorio:

ST/ - Scripts para el análisis de datos de ST (*Visium*) con *Giotto*:

-   00_Construir_objeto.R: Construcción del objeto *Giotto* a partir de los datos descargados.

-   01_Control_calidad.R: Cálculo de métricas de calidad de *spots* y genes.

-   02_Normalizacion.R: Normalización de los conteos.

-   03_HVGs.R: Identificación de genes altamente variables.

-   04_RedDim.R: Reducción de dimensionalidades.

-   05_Leiden.R: Agrupamiento de *spots* mediante el algoritmo de Leiden.

-   06_Deconvolucion.R: Deconvolución utilizandos datos de snRNA-seq.

-   07_Agrupamiento_espacial.R: Agrupamiento espacial de *spots* con modelos HMRF.

-   07_1_Delimitar_regiones.R: Clasificación de dominios espaciales en sustancia blanca y sustancia gris.

-   08_Scran.R: Análisis de expresión diferencial con el método Scran.

snRNA-seq/ - Scripts para el análisis de datos de snRNA-seq:

-   00_Construir_SCE.R: Construcción del objeto *SingleCellExperiment* a partir de los datos crudos.

-   01_Control_calidad.R: Cálculo de métricas de calidad por núcleo.

-   02_Normalizacion.R: Conversión a objeto *Giotto* y normalización de los conteos.

-   03_HVGs_RedDim.R: Selección de genes altamente variables y reducción de dimensionalidad.

Los datos no se incluyen en el respositorio, pero pueden descargarse desde GEO ([GSE233208](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233208)).
