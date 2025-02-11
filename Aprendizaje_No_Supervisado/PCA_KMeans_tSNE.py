##############################################################
## Código de aprendizaje no supervisado por K-means y t-SNE ##
##############################################################

######################################
# Cargar los datos de LFC y p-valores #
######################################

import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# Cargar los datos de LFC y p-valores
lfc_file_path = 'input-values.csv'
#########pval_file_path = 'input-p-value.csv'

lfc_data = pd.read_csv(lfc_file_path, delimiter=';')
#pval_data = pd.read_csv(pval_file_path, delimiter=';')

# Extraer nombres de tejidos y datos de expresión genética
gene_names = lfc_data['Unnamed: 0']  # Nombres de genes
lfc_values = lfc_data.iloc[:, 1:]      # Expresión genética (valores LFC)

# Eliminar datos no numéricos y faltantes
lfc_values = lfc_values.apply(pd.to_numeric, errors='coerce').dropna(axis=1)

# Verificar columnas con varianza cero
zero_variance_cols = lfc_values.var()[lfc_values.var() == 0].index
lfc_values = lfc_values.drop(columns=zero_variance_cols)

# Estandarizar los valores LFC para clustering
scaler = StandardScaler()
lfc_scaled = scaler.fit_transform(lfc_values)

# Determinar el número óptimo de clústeres usando el método del codo
inertia = []
cluster_range = range(1, 8)
for n_clusters in cluster_range:
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    kmeans.fit(lfc_scaled)
    inertia.append(kmeans.inertia_)

# Graficar la curva del codo
plt.figure(figsize=(4, 5))
plt.plot(cluster_range, inertia, marker='o', linestyle='--')
plt.title('Método del codo para los clústeres óptimos')
plt.xlabel('Número de clústeres')
plt.ylabel('Inercia')
plt.show()

##########################################################
# Realizar clustering con el número óptimo de clústeres #
##########################################################

optimal_clusters = 2  # Ajustar según el resultado del método del codo
kmeans = KMeans(n_clusters=optimal_clusters, random_state=30, n_init=10)
clusters = kmeans.fit_predict(lfc_scaled)

# Agregar asignaciones de clúster al DataFrame de LFC
lfc_data['Cluster'] = clusters

# Guardar los datos agrupados para revisión
lfc_data.to_csv('clustered_data.csv', index=False)
print("Datos agrupados guardados en 'clustered_data.csv'")

import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# Reducir dimensiones a 2D usando PCA
pca = PCA(n_components=2)
lfc_pca = pca.fit_transform(lfc_scaled)

# Graficar los resultados del clustering en la proyección PCA
plt.figure(figsize=(8, 5))
for cluster in range(optimal_clusters):
    cluster_points = lfc_pca[clusters == cluster]
    plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f"Clúster {cluster + 1}", s=100, alpha=0.7)

    # Calcular el centroide del clúster
    centroid = np.mean(cluster_points, axis=0)
    # Dibujar un círculo alrededor del centroide del clúster
    circle = plt.Circle(centroid, 1.9, color='black', fill=False, linestyle='--')
    plt.gca().add_patch(circle)

    # Etiquetar cada punto en el clúster con el gen correspondiente
    for i, gene in enumerate(gene_names[clusters == cluster]):
        plt.text(cluster_points[i, 0], cluster_points[i, 1], gene, fontsize=9, ha='right', alpha=0.7)

plt.title("Clústeres (Proyección de PCA)")
plt.xlabel("Componente 1 de PCA (55.59%)")
plt.ylabel("Componente 2 de PCA (21.69%)")
plt.legend()
plt.show()

# Opcional: Usar t-SNE para visualización (mejor para separaciones no lineales)
tsne = TSNE(n_components=2, random_state=1, perplexity=4, n_iter=300)
lfc_tsne = tsne.fit_transform(lfc_scaled)

plt.figure(figsize=(8, 5))
for cluster in range(optimal_clusters):
    cluster_points = lfc_tsne[clusters == cluster]
    plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f"Clúster {cluster + 1}", s=100, alpha=0.7)

    # Calcular el centroide del clúster
    centroid = np.mean(cluster_points, axis=0)
    # Dibujar un círculo alrededor del centroide del clúster
    circle = plt.Circle(centroid, 26.3, color='black', fill=False, linestyle='--')
    plt.gca().add_patch(circle)

    # Etiquetar cada punto en el clúster con el gen correspondiente
    for i, gene in enumerate(gene_names[clusters == cluster]):
        plt.text(cluster_points[i, 0], cluster_points[i, 1], gene, fontsize=9, ha='right', alpha=0.7)

plt.title("Clústeres (Proyección de t-SNE)")
plt.xlabel("Componente 1 de t-SNE")
plt.ylabel("Componente 2 de t-SNE")
plt.legend()
plt.show()

# Graficar los resultados del clustering en la proyección PCA
plt.figure(figsize=(8, 5))
centroids = []

# Graficar puntos y calcular centroides para cada clúster
for cluster in range(optimal_clusters):
    cluster_points = lfc_pca[clusters == cluster]
    plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f"Clúster {cluster + 1}", s=100, alpha=0.7)

    # Calcular y almacenar el centroide del clúster
    centroid = np.mean(cluster_points, axis=0)
    centroids.append(centroid)

    # Etiquetar cada punto en el clúster con el gen correspondiente
    for i, gene in enumerate(gene_names[clusters == cluster]):
        plt.text(cluster_points[i, 0], cluster_points[i, 1], gene, fontsize=9, ha='right', alpha=0.7)

# Extraer centroides para los dos clústeres
centroid1, centroid2 = centroids

# Calcular el punto medio entre los centroides
midpoint = (centroid1 + centroid2) / 2

# Calcular la pendiente de la bisectriz perpendicular
delta_y = centroid2[1] - centroid1[1]
delta_x = centroid2[0] - centroid1[0]
if delta_x != 0:
    slope_perpendicular = -delta_x / delta_y
else:
    slope_perpendicular = float('inf')

# Definir la línea de la frontera de decisión
x_vals = np.array(plt.gca().get_xlim())
if slope_perpendicular != float('inf'):
    y_vals = midpoint[1] + slope_perpendicular * (x_vals - midpoint[0])
else:
    x_vals = np.full(2, midpoint[0])
    y_vals = np.array(plt.gca().get_ylim())

# Graficar la frontera de decisión
plt.plot(x_vals, y_vals, 'k--')

# Establecer límites en el eje y
plt.ylim(-3.5, 3.5)

plt.title("Clústeres (Proyección de PCA)")
plt.xlabel("Componente 1 de PCA (55.59%)")
plt.ylabel("Componente 2 de PCA (21.69%)")
plt.legend()
plt.show()

# Agregar asignaciones de clúster al conjunto de datos original
lfc_data['Cluster'] = clusters

# Guardar los resultados del clustering en un archivo CSV
output_file_path = 'clustered_results.csv'
lfc_data.to_csv(output_file_path, index=False)

print(f"Resultados del clustering guardados en {output_file_path}")

#######################
# Graficar Scree plot #
#######################

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Realizar PCA sobre los datos de LFC
pca = PCA()
lfc_pca = pca.fit_transform(lfc_scaled)

# Extraer proporción de varianza explicada y varianza acumulada
explained_variance_ratio = pca.explained_variance_ratio_ * 100  # Convertir a porcentaje
cumulative_variance = np.cumsum(explained_variance_ratio)

# Imprimir proporción de varianza explicada y varianza acumulada por componente
print("Componente Principal | Varianza Explicada (%) | Varianza Acumulada (%)")
for i, (evr, cv) in enumerate(zip(explained_variance_ratio, cumulative_variance), 1):
    print(f"Componente {i}: {evr:.2f}% | Cumulada: {cv:.2f}%")

# Crear el Scree plot
plt.figure(figsize=(4, 5))
components = range(1, len(explained_variance_ratio) + 1)

# Graficar la varianza explicada como barras azules
plt.bar(components, explained_variance_ratio, color='#06437a', alpha=0.7, label='Varianza explicada (%)')

# Graficar la varianza acumulada como una línea roja
plt.plot(components, cumulative_variance, color='#f67f47', marker='o', linestyle='--', label='Varianza acumulada (%)')

# Agregar etiquetas y título
plt.xlabel('Componente principal')
plt.ylabel('Varianza explicada (%)')
plt.title('Scree Plot de los resultados de PCA')
plt.xticks(components)  # Asegurar que todos los componentes estén etiquetados
plt.legend()
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Mostrar el gráfico
plt.tight_layout()
plt.show()
