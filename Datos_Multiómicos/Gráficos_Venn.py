####################################################
## Análisis de coincidencias por gráficos de Venn ##
####################################################

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from itertools import combinations
import os
import zipfile
import matplotlib.patches as patches
import pandas as pd

# Recargar el archivo CSV para restaurar el contexto
clusters_data = pd.read_csv("Clusters.csv", sep=';')
clusters_data.head()

# Función para calcular el índice de Jaccard
def jaccard_index(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0

# Función auxiliar para guardar y comprimir diagramas de Venn
def save_and_compress_images(image_folder, zip_filename):
    with zipfile.ZipFile(zip_filename, 'w') as zipf:
        for root, _, files in os.walk(image_folder):
            for file in files:
                file_path = os.path.join(root, file)
                zipf.write(file_path, os.path.relpath(file_path, image_folder))

# Función auxiliar para trazar diagramas de Venn para dos conjuntos
def plot_venn2(set1, set2, set_labels, title, image_path):
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    # Crear diagrama base usando elipses
    colors = ['white', 'white']
    sets = [set1, set2]
    centers = [(0.4, 0.5), (0.6, 0.5)]
    width, height = 0.3, 0.2

    for idx, center in enumerate(centers):
        ellipse = patches.Ellipse(center, width, height, color=colors[idx], alpha=0.4, label=set_labels[idx])
        ax.add_patch(ellipse)

    # Anotar los elementos de los conjuntos basados en las intersecciones
    unique_elements = [set1, set2, set1 & set2]
    positions = [(0.3, 0.5), (0.7, 0.5), (0.5, 0.5)]

    for idx, elements in enumerate(unique_elements):
        if elements:
            plt.text(positions[idx][0], positions[idx][1], "\n".join(map(str, elements)), ha='center', va='center', fontsize=8, bbox=dict(facecolor='white', alpha=0.7))

    plt.title(title)
    plt.legend(loc='upper right')
    plt.axis('off')
    plt.savefig(image_path)
    plt.close()

# Función auxiliar para trazar un Venn personalizado para tres conjuntos
def plot_venn3(set1, set2, set3, set_labels, title, image_path):
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    # Crear diagrama base de Venn usando elipses
    colors = ['white', 'white', 'white']
    sets = [set1, set2, set3]
    centers = [(0.4, 0.6), (0.6, 0.6), (0.5, 0.4)]
    width, height = 0.3, 0.2

    for idx, center in enumerate(centers):
        ellipse = patches.Ellipse(center, width, height, color=colors[idx], alpha=0.4, label=set_labels[idx])
        ax.add_patch(ellipse)

    # Anotar los elementos de los conjuntos basados en las intersecciones
    unique_elements = [set1, set2, set3, set1 & set2, set1 & set3, set2 & set3, set1 & set2 & set3]
    positions = [(0.3, 0.6), (0.7, 0.6), (0.5, 0.3), (0.5, 0.7), (0.4, 0.5), (0.6, 0.5), (0.5, 0.5)]

    for idx, elements in enumerate(unique_elements):
        if elements:
            plt.text(positions[idx][0], positions[idx][1], "\n".join(map(str, elements)), ha='center', va='center', fontsize=8, bbox=dict(facecolor='white', alpha=0.7))

    plt.title(title)
    plt.legend(loc='upper right')
    plt.axis('off')
    plt.savefig(image_path)
    plt.close()

# Función auxiliar para trazar un Venn personalizado para cuatro conjuntos
def plot_venn4(set1, set2, set3, set4, set_labels, title, image_path):
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    # Crear diagrama base de Venn usando elipses
    colors = ['white', 'white', 'white', 'white']
    sets = [set1, set2, set3, set4]
    centers = [(0.3, 0.5), (0.5, 0.5), (0.4, 0.7), (0.4, 0.3)]
    width, height = 0.3, 0.2

    for idx, center in enumerate(centers):
        ellipse = patches.Ellipse(center, width, height, color=colors[idx], alpha=0.4, label=set_labels[idx])
        ax.add_patch(ellipse)

    # Anotar los elementos de los conjuntos basados en las intersecciones
    all_intersections = [set1 & set2, set1 & set3, set1 & set4, set2 & set3, set2 & set4, set3 & set4, set1 & set2 & set3, set1 & set2 & set4, set1 & set3 & set4, set2 & set3 & set4, set1 & set2 & set3 & set4]
    unique_elements = [set1, set2, set3, set4] + all_intersections

    for idx, elements in enumerate(unique_elements):
        if elements:
            pos_x = centers[idx % 4][0] if idx < 4 else 0.5
            pos_y = centers[idx % 4][1] if idx < 4 else 0.5
            plt.text(pos_x, pos_y, "\n".join(map(str, elements)), ha='center', va='center', fontsize=8, bbox=dict(facecolor='white', alpha=0.7))

    plt.title(title)
    plt.legend(loc='upper right')
    plt.axis('off')

    plt.savefig(image_path)
    plt.close()

# Función auxiliar para procesar los clusters y trazar diagramas de Venn
def process_clusters(cluster_data, columns, title, venn_type='venn2', output_folder='venn_outputs'):
    # Crear la carpeta de salida si no existe
    os.makedirs(output_folder, exist_ok=True)

    # Extraer los números de los clusters y agrupar los genes por clusters para cada columna
    cluster_sets = {
        col: {
            cluster: set(cluster_data.loc[cluster_data[col] == cluster, 'Genes'])
            for cluster in cluster_data[col].dropna().unique()
        }
        for col in columns
    }

    # Generar diagramas de Venn y calcular el índice de Jaccard para cada número de cluster
    results = {}
    for cluster_num in set().union(*[set(cluster_sets[col].keys()) for col in columns]):
        # Obtener los conjuntos de genes para este cluster a través de las columnas
        sets = [cluster_sets[col].get(cluster_num, set()) for col in columns]

        # Omitir clusters sin solapamientos o valores faltantes
        if venn_type == 'venn2' and len(sets) == 2:
            set1, set2 = sets
            image_path = os.path.join(output_folder, f"{title.replace(' ', '_')}_Cluster_{cluster_num}.png")
            plot_venn2(set1, set2, columns, f"{title} - Cluster {cluster_num}", image_path)
            results[cluster_num] = jaccard_index(set1, set2)

        elif venn_type == 'venn3' and len(sets) == 3:
            set1, set2, set3 = sets
            image_path = os.path.join(output_folder, f"{title.replace(' ', '_')}_Cluster_{cluster_num}.png")
            plot_venn3(set1, set2, set3, columns, f"{title} - Cluster {cluster_num}", image_path)
            results[cluster_num] = {
                f"{columns[i]} & {columns[j]}": jaccard_index(sets[i], sets[j])
                for i, j in combinations(range(3), 2)
            }

        elif venn_type == 'venn4' and len(sets) == 4:
            set1, set2, set3, set4 = sets
            image_path = os.path.join(output_folder, f"{title.replace(' ', '_')}_Cluster_{cluster_num}.png")
            plot_venn4(set1, set2, set3, set4, columns, f"{title} - Cluster {cluster_num}", image_path)
            # Calcular los índices de Jaccard para todas las combinaciones por pares
            results[cluster_num] = {
                f"{columns[i]} & {columns[j]}": jaccard_index(sets[i], sets[j])
                for i, j in combinations(range(4), 2)
            }

    return results

# Comprimir las imágenes de salida en un archivo ZIP
def compress_outputs(output_folder='venn_outputs', zip_filename='Venn_Diagrams.zip'):
    save_and_compress_images(output_folder, zip_filename)
    return zip_filename

############################
# Generar Gráficos de Venn #
############################

# Tarea 1: Transcriptómica del cáncer (columnas 3 y 6)
expr_fisiologica_results = process_clusters(
    clusters_data,
    ['Transcriptómica - Fisiológico', 'Transcriptómica - Cáncer'],
    title="Transcriptómica del cáncer",
    venn_type='venn2'
)

expr_fisiologica_results

# Tarea 2: Proteómica del cáncer (columnas 3 y 6)
expr_fisiologica_results = process_clusters(
    clusters_data,
    ['Proteómica - Fisiológico', 'Proteómica - Cáncer'],
    title="Proteómica del cáncer",
    venn_type='venn2'
)

expr_fisiologica_results

# Tarea 3: Expresión tumoral multiómica (columnas 2, 4 y 5)
expr_tumoral_results = process_clusters(
    clusters_data,
    ['Genómica - Cáncer', 'Transcriptómica - Cáncer', 'Proteómica - Cáncer'],
    title="Expresión tumoral multiómica",
    venn_type='venn3'
)

expr_tumoral_results


# Tarea 4: Expresión fisiológica (columnas 3 y 4)
cambio_expr_tumores_results = process_clusters(
    clusters_data,
    ['Transcriptómica - Fisiológico', 'Proteómica - Fisiológico'],
    title="Expresión fisiológica",
    venn_type='venn2'
)

cambio_expr_tumores_results

# Tarea 5: Expresión diferencial multiómica (columnas 3 y 4)
cambio_expr_omicas_results = process_clusters(
    clusters_data,
    ['Transcriptómica - Fisiológico', 'Transcriptómica - Cáncer', 'Proteómica - Fisiológico', 'Proteómica - Cáncer'],
    title="Expresión diferencial multiómica",
    venn_type='venn4'
)

cambio_expr_omicas_results

# Especificar la carpeta de salida donde se guardaron las imágenes de los diagramas de Venn
output_folder = 'venn_outputs'

# Especificar el nombre del archivo ZIP a crear
zip_filename = 'Venn_Diagrams.zip'

# Llamar a la función compress_outputs para crear el archivo ZIP
zip_file_path = compress_outputs(output_folder=output_folder, zip_filename=zip_filename)

# Imprimir confirmación de la creación del archivo ZIP
print(f"Archivo ZIP creado: {zip_file_path}")
