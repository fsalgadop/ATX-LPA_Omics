
#####################################################
## Correlación entre las bases de datos de RNA-seq ##
#####################################################

#############################
## Cargar las bibliotecas requeridas ##
#############################
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ggplot2)
library(pheatmap)
library(reshape2)
library(dplyr)

##############################
## Cargar Datos y Preprocesar ##
##############################

# Establecer el directorio que contiene los archivos CSV
data_dir <- "" # Añadir el directorio con los datos de RNA-seq de uno o todos los tejidos estudiados (P.ej. "Pulmón")

# Definir los genes de interés
genes_of_interest <- c("ENPP2","ITGA2B","ITGA5","ITGA8","ITGAV","ITGB1","ITGB3",
                       "ITGB5","ITGB6","ITGB8","LPAR1","LPAR2","LPAR3","LPAR4",
                       "LPAR5","LPAR6","PLPP1","PLPP2","PLPP3")

# Listar todos los archivos CSV en el directorio
file_list <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

# Inicializar una lista vacía para almacenar los dataframes
data_list <- list()

# Leer todos los archivos CSV y almacenarlos en una lista
for (file in file_list) {
  df <- read.csv(file, row.names = 1, sep = ";")  # Suponiendo que la primera columna contiene los nombres de las filas (genes)
  # Filtrar para incluir solo los genes de interés
  df_filtered <- df[genes_of_interest, "log2FoldChange", drop = FALSE] # Mantener solo la columna log2FoldChange
  data_list[[file]] <- df_filtered
}

# Asegurar la alineación de filas ordenando los nombres de las filas de manera consistente entre los datasets
common_genes <- Reduce(intersect, lapply(data_list, rownames))  # Encontrar los genes comunes en todos los datasets
# Alinear los datasets por los genes comunes
aligned_data <- lapply(data_list, function(df) df[common_genes, , drop = FALSE])
# Combinar los datasets por columnas (suponiendo que la columna "log2FoldChange" es la única de interés)
merged_data <- do.call(cbind, aligned_data)
# Etiquetar las columnas para identificar los datasets de origen
colnames(merged_data) <- paste0("Dataset_", seq_along(aligned_data))
# Valores de Log2FoldChange para cada gen en los diferentes datasets
log2FoldChange_data <- as.matrix(merged_data)
# Verificar la estructura de los datos finales
head(log2FoldChange_data)

####################################
## Calcular Correlaciones de Pearson ##
####################################

log2FoldChange_data <- as.data.frame(log2FoldChange_data)

# Función para eliminar filas con valores NA de un dataframe
remove_na_rows <- function(df) {
  df_clean <- df[complete.cases(df), ]  # Mantener solo las filas sin valores NA
  return(df_clean)
}

log2FoldChange_data <- remove_na_rows(log2FoldChange_data)

# Calcular la matriz de correlación de Pearson
cor_matrix <- cor(t(log2FoldChange_data), method = "pearson")

# Ordenar las filas y columnas alfabéticamente
gene_order <- order(rownames(cor_matrix))
sorted_cor_matrix <- cor_matrix[gene_order, gene_order]

# Convertir la matriz de R2 ordenada en un dataframe para facilitar su visualización
cor_table <- as.data.frame(sorted_cor_matrix)

# Mostrar la tabla R2
print(cor_table)

# Graficar el Heatmap de Correlaciones
pheatmap(cor_matrix,
         display_numbers = TRUE,  # Mostrar los valores de la correlación de Pearson dentro de las celdas
         number_format = "%.2f",  # Formatear los números con 2 decimales
         main = "Matriz de correlación - Pulmón")

##################################
## Gráfico de dispersión de Correlaciones ##
##################################

# Gráfico de dispersión de ENPP2 vs otros genes en el primer dataset (ejemplo de gráficos de dispersión)
enpp2_vs_others <- log2FoldChange_data["ENPP2", ]
df_scatter <- data.frame(ENPP2 = enpp2_vs_others)

# Añadir otros genes al dataframe para graficar
for (gene in genes_of_interest[genes_of_interest != "ENPP2"]) {
  df_scatter[[gene]] <- log2FoldChange_data[gene, ]
}

# Convertir a formato largo para ggplot
df_long <- melt(df_scatter, id.vars = "ENPP2")

# Gráfico de dispersión con ggplot
ggplot(df_long, aes(x = ENPP2, y = value, color = variable)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "ENPP2 log2FoldChange", y = "Otros Genes log2FoldChange",
       title = "Gráfico de Dispersión: ENPP2 vs Otros Genes") +
  theme_minimal()

####################################
## Correlaciones Específicas de ENPP2 ##
####################################

# 1. Correlaciones de ENPP2 vs otros genes:
# Verificar qué genes de genes_of_interest están presentes en log2FoldChange_data
missing_genes <- setdiff(genes_of_interest, rownames(log2FoldChange_data))

if (length(missing_genes) > 0) {
  print("Los siguientes genes están ausentes en los datos:")
  print(missing_genes)
} else {
  print("Todos los genes están presentes.")
}

# Solo proceder con los genes que están presentes en los datos
present_genes <- intersect(genes_of_interest, rownames(log2FoldChange_data))

# Calcular la correlación de Pearson de ENPP2 con cada uno de los otros genes presentes
cor_enpp2_vs_others <- sapply(present_genes[present_genes != "ENPP2"], function(gene) {
  cor(log2FoldChange_data["ENPP2", ], log2FoldChange_data[gene, ], use = "pairwise.complete.obs")
})

# Mostrar las correlaciones
cor_enpp2_vs_others


# 2. Correlación de ENPP2 y otros genes individuales:

# Definir grupos
group_LPA1_3 <- c("LPAR1", "LPAR2", "LPAR3")
group_LPA4_6 <- c("LPAR4", "LPAR5", "LPAR6")
group_LPPs   <- c("PLPP1", "PLPP2", "PLPP3")

# Combinar todos los genes de los grupos en una lista para verificar su presencia
all_group_genes <- c(group_LPA1_3, group_LPA4_6, group_LPPs)

# Verificar si faltan genes en los datos
missing_group_genes <- setdiff(all_group_genes, rownames(log2FoldChange_data))

if (length(missing_group_genes) > 0) {
  print("Los siguientes genes de los grupos están ausentes en los datos:")
  print(missing_group_genes)
} else {
  print("Todos los genes de los grupos están presentes.")
}

# Mantener solo los genes que están presentes en los datos
present_group_LPA1_3 <- intersect(group_LPA1_3, rownames(log2FoldChange_data))
present_group_LPA4_6 <- intersect(group_LPA4_6, rownames(log2FoldChange_data))
present_group_LPPs   <- intersect(group_LPPs, rownames(log2FoldChange_data))

# Calcular la correlación entre ENPP2 y cada grupo de genes
cor_enpp2_vs_LPA1_3 <- sapply(present_group_LPA1_3, function(gene) {
  cor(log2FoldChange_data["ENPP2", ], log2FoldChange_data[gene, ], use = "pairwise.complete.obs")
})

cor_enpp2_vs_LPA4_6 <- sapply(present_group_LPA4_6, function(gene) {
  cor(log2FoldChange_data["ENPP2", ], log2FoldChange_data[gene, ], use = "pairwise.complete.obs")
})

cor_enpp2_vs_LPPs <- sapply(present_group_LPPs, function(gene) {
  cor(log2FoldChange_data["ENPP2", ], log2FoldChange_data[gene, ], use = "pairwise.complete.obs")
})

# Mostrar las correlaciones para cada grupo
list(
  ENPP2_vs_LPA1_3 = cor_enpp2_vs_LPA1_3,
  ENPP2_vs_LPA4_6 = cor_enpp2_vs_LPA4_6,
  ENPP2_vs_LPPs = cor_enpp2_vs_LPPs
)


# 3. Regresión lineal

# Función para realizar regresión lineal para cada gen
plot_regression_for_genes <- function(target_gene, all_genes) {
  # Verificar si el gen objetivo está en los datos
  if (!(target_gene %in% rownames(log2FoldChange_data))) {
    warning(paste("El gen objetivo", target_gene, "no se encuentra en los datos."))
    return(NULL)
  }
  
  # Crear una lista vacía para almacenar los gráficos
  plots <- list()
  
  # Recorrer todos los genes (excluyendo el gen objetivo)
  for (gene in all_genes) {
    if (gene != target_gene && gene %in% rownames(log2FoldChange_data)) {
      # Preparar un dataframe para la regresión lineal
      data_for_regression <- data.frame(
        TargetGene = log2FoldChange_data[target_gene, ],
        Gene = log2FoldChange_data[gene, ]
      )
      
      # Eliminar filas con NA
      data_for_regression <- na.omit(data_for_regression)
      
      # Verificar si tenemos suficientes datos para la regresión
      if (nrow(data_for_regression) < 2) {
        warning(paste("No hay suficientes datos para la regresión entre", target_gene, "y", gene))
        next
      }
      
      # Ajustar un modelo de regresión lineal
      model <- lm(Gene ~ TargetGene, data = data_for_regression)
      
      # Imprimir el resumen del modelo de regresión lineal
      print
