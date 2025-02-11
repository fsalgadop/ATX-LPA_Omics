##############################################################
## Gráfico de análisis de supervivencia con datos genómicos ##
##############################################################

import matplotlib.pyplot as plt
import pandas as pd

# Datos genómicos
data = {
    'Gene': ['PLPP3', 'PLPP2', 'PLPP1', 'LPAR6', 'LPAR5', 'LPAR4', 'LPAR3', 'LPAR2', 'LPAR1',
             'ITGB8', 'ITGB6', 'ITGB5', 'ITGB3', 'ITGB1', 'ITGAV', 'ITGA8', 'ITGA5', 'ITGA2B',
             'ENPP2', 'Alteraciones'],
    'HR': [1.219, 1.045, 0.915, 0.974, 1.42, 0.95, 0.576, 1.048, 1.381, 1.252, 1.21, 1.177, 1.073,
           1.369, 1.065, 1.217, 0.944, 1.177, 1.027, 1.115],
    'Error': [0.721, 0.535, 0.535, 0.408, 0.839, 0.593, 0.392, 0.728, 0.994, 0.65, 0.797, 0.648,
              0.736, 1.114, 0.553, 0.555, 0.615, 0.847, 0.296, 0.162],
    'p': [0.101, 0.559, 0.628, 0.784, 0.00184, 0.226, 0.0117, 0.33, 0.297, 0.0118, 0.367, 0.204,
          0.27, 0.013, 0.761, 0.133, 0.422, 0.294, 0.826, 0.002703]
}


# Crear DataFrame
df = pd.DataFrame(data)

# Gráfica
plt.figure(figsize=(6, 7))

# Graficar los puntos con barras de error
for i, row in df.iterrows():
    # Establecer color en rojo si el valor p < 0.05, de lo contrario negro
    point_color = 'red' if row['p'] < 0.05 else 'black'

    plt.errorbar(row['HR'], i, xerr=row['Error'], fmt='o', color=point_color,
                 ecolor='black', elinewidth=1, capsize=3)  # Barras de error en negro y más delgadas

    # Agregar valores p a la derecha
    plt.text(row['HR'] + row['Error'] + 0.05, i, f"p={row['p']:.4f}", va='center')

# Dibujar una línea vertical discontinua en X=1
plt.axvline(x=1, color='black', linestyle='dashed', linewidth=1)

# Etiquetas y título
plt.title("Análisis de supervivencia")
plt.xlabel("Factor de riesgo")
plt.ylabel("Genes")

# Establecer límites del eje X
plt.xlim(0, 3)

# Establecer las marcas del eje Y con los nombres ordenados de los genes
plt.yticks(range(len(df)), df['Gene'])

plt.grid(True)
plt.show()
