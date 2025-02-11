#####################################################
## Análisis de LFC de matrices de datos de RNA-seq ##
#####################################################

#############################
## Cargar bibliotecas requeridas ##
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
## Cargar datos y preprocesar ##
##############################

# Definir genes de interés
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

# Ordenar las filas y columnas alfabéticamente
gene_order <- order(rownames(cor_matrix))
sorted_cor_matrix <- cor_matrix[gene_order, gene_order]

# Convertir la matriz R2 ordenada a un data frame para facilitar la visualización
cor_table <- as.data.frame(sorted_cor_matrix)

# Mostrar la tabla R2
print(cor_table)

# Definir la paleta de colores personalizada
custom_colors <- colorRampPalette(c("white", "#2E6F40"))(120)

# Generar el mapa de calor con el esquema de colores personalizado
pheatmap(cor_matrix, 
         main = "Expresión diferencial", 
         color = custom_colors, 
         breaks = seq(-5, +10, length.out = 101))  # Ajustar los puntos de corte para el rango -10 a +5

custom_colors <- colorRampPalette(c("lightskyblue3", "white", "tomato3"))(100)

pheatmap(cor_matrix, 
         main = "Expresión diferencial", 
         color = custom_colors, 
         breaks = seq(-6, +6, length.out = 100))  # Ajustar los puntos de corte para el rango -10 a +5

# Importar la matriz de p-valores
p_value_matrix <- read.delim("C:/Users/fodag/OneDrive/Escritorio/UNIR/Master en Bioinformática/TFM/Databases/To process/Example tissues_old/LFC_LPPs - p-value - stars.csv", sep = ";")
rownames(p_value_matrix) <- p_value_matrix$X
p_value_matrix <- p_value_matrix[,-1]
p_value_matrix <- as.matrix(p_value_matrix)

# Crear el mapa de calor y añadir los p-valores como texto en las celdas
pheatmap(cor_matrix, 
         display_numbers = format(p_value_matrix, digits = 2),  # Añadir los p-valores como texto
         number_format = "%.1f",  # Formatear los p-valores con 2 decimales
         main = "Expresión diferencial",
         color = custom_colors,
         breaks = seq(-6, +6, length.out = 100)  # Ajustar los puntos de corte para el rango -10 a +5
         )
