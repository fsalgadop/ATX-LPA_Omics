####################################
## Análisis de expresión proteica ##
####################################

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
## Cargar y preprocesar los datos ##
##############################

# Definir los genes de interés
genes_of_interest <- c("ENPP2","ITGA2B","ITGA5","ITGA8","ITGAV","ITGB1","ITGB3",
                       "ITGB5","ITGB6","ITGB8","LPAR1","LPAR2","LPAR3","LPAR4",
                       "LPAR5","LPAR6","PLPP1","PLPP2","PLPP3")

cor_matrix <- read.delim("", sep = ";")
rownames(cor_matrix) <- cor_matrix$X
cor_matrix <- cor_matrix[,-1]
cor_matrix <- as.matrix(cor_matrix)

######################
## Crear mapa de calor ##
######################

# Ordenar las filas y las columnas alfabéticamente
gene_order <- order(rownames(cor_matrix))
sorted_cor_matrix <- cor_matrix[gene_order, gene_order]

# Convertir la matriz R2 ordenada en un data frame para una visualización más fácil
cor_table <- as.data.frame(sorted_cor_matrix)

# Mostrar la tabla R2
print(cor_table)

# Definir la paleta de colores personalizada
custom_colors <- colorRampPalette(c("white", "#2E6F40"))(100)

# Generar el mapa de calor con el esquema de colores personalizado
custom_colors <- colorRampPalette(c("white", "#2E6F40"))(100)
pheatmap(cor_matrix, 
         main = "Niveles de expresión de proteínas", 
         color = custom_colors, 
         breaks = seq(0, +3, length.out = 100))  # Ajustar los puntos de interrupción para el rango de -10 a +5

custom_colors <- colorRampPalette(c("white", "#2E6F40"))(100)
pheatmap(cor_matrix, 
         main = "Niveles de expresión en cáncer - HPA", 
         color = custom_colors, 
         breaks = seq(0, +3, length.out = 100))  # Ajustar los puntos de interrupción para el rango de -10 a +5
