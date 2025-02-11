##########################################
## Gráfico_Alteraciones_Genéticas_ENPP2 ##
##########################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Datos de las alteraciones genómicas de ENPP2 (sin normalización al 100%)
data = {
    "Tumor": ["Cáncer de vejiga", "Cáncer de mama", "Cáncer colorrectal", "Cáncer endometrial",
              "Cáncer esofagogástrico", "Cáncer hepatobiliar", "Melanoma",
              "Cáncer de pulmón no microcítico", "Tumor epitelial ovárico", "Cáncer pancreático",
              "Cáncer de próstata"],
    "Amplificación": [5.109489051, 11.71586716, 3.03030303, 4.266211604, 8.360128617, 10.21505376,
                      3.153153153, 3.133903134, 21.23287671, 8.152173913, 6.882591093],
    "Delección": [0, 0.184501845, 0, 0, 0, 0.268817204, 0, 0, 0.51369863, 0, 0.20242915],
    "Múltiples alteraciones": [0, 0.36900369, 0.168350168, 0, 0, 0.537634409, 0, 0.094966762, 0.51369863,
                               0, 0],
    "Mutación": [1.946472019, 0.922509225, 2.525252525, 5.972696246, 1.92926045, 0, 4.279279279,
                 3.51377018, 0.856164384, 0.543478261, 0.4048583]
}

# Crear un DataFrame
df = pd.DataFrame(data)

# Graficar
fig, ax = plt.subplots(figsize=(3, 4))

# Colores para cada tipo de alteración
colors = {
    "Amplificación": "coral",
    "Delección": "darkcyan",
    "Múltiples alteraciones": "slategray",
    "Mutación": "darkseagreen"
}

# Gráfico de barras apiladas con datos reales (sin normalización)
ax.barh(df["Tumor"], df["Amplificación"], label="Amplificación", color=colors["Amplificación"])
ax.barh(df["Tumor"], df["Delección"], left=df["Amplificación"], label="Delección", color=colors["Delección"])
ax.barh(df["Tumor"], df["Múltiples alteraciones"], left=df["Amplificación"] + df["Delección"],
        label="Múltiples alteraciones", color=colors["Múltiples alteraciones"])
ax.barh(df["Tumor"], df["Mutación"], left=df["Amplificación"] + df["Delección"] + df["Múltiples alteraciones"],
        label="Mutación", color=colors["Mutación"])

# Etiquetas y personalización
ax.set_xlabel("Frecuencia (%)")
ax.set_title("Alteraciones de ENPP2 en tumores")
ax.legend(title="Tipos de alteraciones", bbox_to_anchor=(1, 1), loc='lower left')
plt.tight_layout()

# Mostrar gráfica
plt.show()
