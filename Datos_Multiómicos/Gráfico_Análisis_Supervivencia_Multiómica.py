################################################################
## Gráfico de análisis de supervivencia con datos multiómicos ##
################################################################

import matplotlib.pyplot as plt
import pandas as pd

# Datos multiómicos
data = {
    'Gene': ['Sin co-ocurrencia tumoral', 'Sin co-ocurrencia fisiológica', 'Co-ocurrencia tumoral',
             'Co-ocurrencia fisiológica', 'Alteraciones'],
    'HR': [1.186, 1.22, 1.162, 1.211, 1.115],
    'Error': [0.101, 0.151, 0.089, 0.1, 0.162],
    'p': [0.145, 0.592, 0.157, 0.0374, 0.002703]
}

# Crear DataFrame
df = pd.DataFrame(data)

# Gráfica
plt.figure(figsize=(6, 3))

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
plt.xlim(0.4, 1.8)

# Establecer las marcas del eje Y con los nombres ordenados de los genes
plt.yticks(range(len(df)), df['Gene'])

plt.grid(True)
plt.show()
