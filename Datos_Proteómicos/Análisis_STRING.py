##############################
## Análisis de Redes STRING ##
##############################

import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import io

# Definir clústeres y sus colores
clusters_genes = {
    0: ['ENPP2','LPAR1','LPAR6','ITGA2B','ITGA5','ITGA8','ITGAV','ITGB1','ITGB3','ITGB5','PLPP1','PLPP3'],
    1: ['LPAR2','LPAR3','LPAR4','LPAR5','ITGB6','ITGB8','PLPP2'],
    2: [''],
    3: ['']
}
cluster_colors = {0: 'skyblue', 1: 'orange', 2: 'lightgreen', 3:'red'}

# Aplanar la lista de genes y crear un mapeo de genes a sus clústeres
gene_to_cluster = {gene: cluster for cluster, genes in clusters_genes.items() for gene in genes}

# Función para consultar la API de STRING y realizar el análisis
def string_analysis(genes, species=9606):
    """
    Consultar la API de STRING para una lista de genes y realizar un análisis de interacción proteína-proteína.

    Parámetros:
    - genes: Lista de nombres de genes.
    - species: ID de taxonomía de NCBI (por defecto 9606 para humanos).

    Retorna:
    - DataFrame de interacciones de STRING.
    """
    string_api_url = "https://string-db.org/api/tsv/network"
    params = {
        "identifiers": "%0d".join(genes),
        "species": species,
        "caller_identity": "your_app_name",  # Reemplazar con tu identificador único de aplicación
    }

    response = requests.get(string_api_url, params=params)
    if response.status_code == 200:
        data = response.text
        string_results = pd.read_csv(io.StringIO(data), sep='\t')
        return string_results
    else:
        print(f"Error con la solicitud a la API de STRING: {response.status_code}")
        return None

# Función para obtener colores de nodos basados en el clúster
def get_node_colors(genes, clusters, cluster_colors):
    """
    Obtener el color de cada nodo basado en su clúster.

    Parámetros:
    - genes: Lista de genes (nodos en la red).
    - clusters: Diccionario que mapea genes a sus clústeres.
    - cluster_colors: Diccionario que mapea IDs de clústeres a colores.

    Retorna:
    - Lista de colores para cada nodo.
    """
    return [cluster_colors[clusters[gene]] if gene in clusters else 'gray' for gene in genes]

# Función para graficar la red de STRING con colores basados en clústeres
def plot_string_network_with_colors(string_data, clusters, cluster_colors, title="STRING Network"):
    """
    Graficar la red de STRING con colores de nodos basados en clústeres.

    Parámetros:
    - string_data: DataFrame de interacciones de STRING.
    - clusters: Diccionario que mapea genes a sus clústeres.
    - cluster_colors: Diccionario que mapea IDs de clústeres a colores.
    - title: Título de la gráfica.
    """
    # Extraer aristas
    edges = string_data[['preferredName_A', 'preferredName_B', 'score']]

    # Crear un grafo de NetworkX
    G = nx.Graph()
    for _, row in edges.iterrows():
        G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])

    # Obtener colores de nodos
    nodes = list(G.nodes)
    node_colors = get_node_colors(nodes, clusters, cluster_colors)

    # Dibujar la red
    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(G, seed=9)  # Distribución de resortes para posicionamiento
    nx.draw_networkx_nodes(G, pos, node_size=150, node_color=node_colors, alpha=1.0)
    nx.draw_networkx_edges(G, pos, width=0.5, alpha=1.0)
    nx.draw_networkx_labels(G, pos, font_size=8, font_color='black')

    # Agregar título y mostrar gráfica
    plt.title(title)
    plt.axis('off')
    plt.tight_layout()
    plt.show()

# Realizar análisis de STRING en todos los genes
all_genes = sum(clusters_genes.values(), [])  # Aplanar las listas de genes
string_results = string_analysis(all_genes)

# Graficar la red de STRING con colores de nodos si los resultados están disponibles
if string_results is not None:
    plot_string_network_with_colors(
        string_results,
        clusters=gene_to_cluster,
        cluster_colors=cluster_colors,
        title="Red STRING: Genes usados para PCA y t-SNE (Coloreados por Clúster)"
    )
