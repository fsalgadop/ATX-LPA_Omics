##########################################
## Gráfico_Alteraciones_Genéticas_Genes ##
##########################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Nuevos datos (manteniendo el orden original)
data = {
    "Gen": ["PLPP3", "PLPP2", "PLPP1", "LPAR6", "LPAR5", "LPAR4", "LPAR3", "LPAR2", "LPAR1", "ITGB8", "ITGB6", "ITGB5",
            "ITGB3", "ITGB1", "ITGAV", "ITGA8", "ITGA5", "ITGA2B", "ENPP2"],
    "Amplificación": [14.88487467, 22.96931623, 9.235052229, 3.768400709, 42.34187185, 5.864630427, 13.73166598, 18.29650994,
                      11.32120548, 27.61908665, 10.69552109, 22.04033852, 13.83772455, 10.78048484, 12.33208481, 16.97407573,
                      20.65696895, 12.15740409, 111.3832313],
    "Deleción": [7.312091304, 20.54049369, 12.36216788, 66.46008386, 2.314915476, 5.332790209, 10.36355038, 4.28120914,
                  3.391869452, 2.750307724, 8.049216791, 5.076702011, 3.361349585, 3.294151914, 16.20370979, 8.71663583,
                  1.008314972, 8.735146956, 3.695342705],
    "Múltiples alteraciones": [0, 1.342210101, 1.503759398, 0.560507031, 0.171232877, 0.194552529, 0.094966762, 0.094966762,
                               0, 0.320191987, 0.168350168, 0.509508641, 0.184501845, 0.191204589, 0.800632068, 0.537576057,
                               0.365785406, 0, 1.683653659],
    "Mutación": [12.22319907, 6.791899686, 9.676370064, 9.111339909, 9.481082788, 21.25072839, 14.79498329, 11.64450934,
                 18.22591211, 29.93211639, 26.80971441, 22.93250788, 32.0241691, 19.62941982, 34.18662965, 52.74727623,
                 28.42361728, 27.27397379, 33.38484496]
}

# Crear un DataFrame
df = pd.DataFrame(data)

# Normalizar los datos al 100% para cada gen
df_normalized = df.drop(columns="Gen").div(df.drop(columns="Gen").sum(axis=1), axis=0) * 100

# Graficar
fig, ax = plt.subplots(figsize=(8, 7))

# Colores para cada tipo de alteración
colors = {
    "Amplificación": "coral",
    "Deleción": "darkcyan",
    "Múltiples alteraciones": "slategray",
    "Mutación": "darkseagreen"
}

# Gráfico de barras apiladas con datos normalizados y ejes invertidos
ax.barh(df["Gen"], df_normalized["Amplificación"], label="Amplificación", color=colors["Amplificación"])
ax.barh(df["Gen"], df_normalized["Deleción"], left=df_normalized["Amplificación"], label="Deleción", color=colors["Deleción"])
ax.barh(df["Gen"], df_normalized["Múltiples alteraciones"], left=df_normalized["Amplificación"] + df_normalized["Deleción"],
        label="Múltiples alteraciones", color=colors["Múltiples alteraciones"])
ax.barh(df["Gen"], df_normalized["Mutación"], left=df_normalized["Amplificación"] + df_normalized["Deleción"] + df_normalized["Múltiples alteraciones"],
        label="Mutación", color=colors["Mutación"])

# Etiquetas y personalización
ax.set_ylabel("Genes")
ax.set_xlabel("Frecuencia (%)")
ax.set_title("Alteraciones genéticas")
ax.legend(title="Tipo de alteraciones", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Mostrar gráfica
plt.show()
