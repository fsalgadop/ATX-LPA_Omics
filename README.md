# ATX-LPA_Omics
Análisis multiómico de la asociación oncológica del eje de señalización ATX-LPA

El código se organiza en las cuatro secciones de resultados, además de una sección de aprendizaje no supervisado.
En cada sección, se incluyen los datos empleados para preparar las figuras del trabajo
Debido a su gran peso, no se adjuntan las bases de datos empleadas de manera directa, sino que se referencian en el trabajo escrito.

# Estructura del Proyecto
**Datos_Transcriptómicos/**
  - Análisis_Correlación_RNA-seq.R
  - RNA-seq_Análisis_LFC.R
  - Expresión_Fisiológica.R
  **Datos/**
    **Valores de expresión fisiológica/**
      - Expresión transcriptómica - HPA_Expresión_Fisiológica.csv
    **Correlación/**
      - Valores_LFC_para_Correlación - Todos los genes.csv
    **LFC_RNA-seq/**
      **Liver/**
        - Liver_2.csv
        - Liver_3.csv
        - Liver_1.csv
      **Lung/**
        - Lung_1.csv
        - Lung_3.csv
        - Lung_2.csv
        - Lung_4.csv
      **Kidney/**
        - Kidney_1.csv
        - Kidney_3.csv
        - Kidney_2.csv
      **Melanoma/**
        - Skin_3.csv
        - Skin_1.csv
        - Skin_2.csv
      **Colon/**
        - Colon_1.csv
        - Colon_3.csv
        - Colon_2.csv
      **Gastric/**
        - Stomach_3.csv
        - Stomach_2.csv
        - Stomach_1.csv
      **Pancreas/**
        - Pancreas_2.csv
        - Pancreas_3.csv
        - Pancreas_1.csv
        - Pancreas_4.csv
      **Breast/**
        - Breast_3.csv
        - Breast_2.csv
        - Breast_1.csv
- README.md
**Datos_Multiómicos/**
  - Gráficos_Venn.py
  **Datos/**
    - Clusters_para_análisis.csv
  - Gráfico_Análisis_Supervivencia_Multiómica.py
**.git/**
  **branches/**
  - config
  **info/**
    - exclude
  - HEAD
  - description
  **logs/**
    - HEAD
    **refs/**
      **remotes/**
        **origin/**
          - HEAD
      **heads/**
        - main
  **refs/**
    **remotes/**
      **origin/**
        - HEAD
    **heads/**
      - main
    **tags/**
  **objects/**
    **info/**
    **pack/**
      - pack-a03cbac5ffd483af6a9ba94189cfb0f87725e0ef.pack
      - pack-a03cbac5ffd483af6a9ba94189cfb0f87725e0ef.idx
  - index
  - packed-refs
  **hooks/**
    - pre-push.sample
    - applypatch-msg.sample
    - commit-msg.sample
    - prepare-commit-msg.sample
    - pre-merge-commit.sample
    - post-update.sample
    - fsmonitor-watchman.sample
    - push-to-checkout.sample
    - update.sample
    - pre-receive.sample
    - pre-commit.sample
    - pre-applypatch.sample
    - pre-rebase.sample
**Aprendizaje_No_Supervisado/**
  - PCA_KMeans_tSNE.py
  **Resultados_Clusterización/**
    **Clusterización - Proteómica_fisiológica/**
      - Elbow.png
      - PCA.png
      **Input/**
        - input-values.csv
      - Scree.png
      - STRING.png
      - clustered_results.csv
    **Clusterización - Genómica/**
      - Elbow.png
      - PCA.png
      **Input/**
        - input-values.csv
      - Scree.png
      - STRING.png
      - clustered_results.csv
    **Clusterización - Transcriptómica_cáncer_RNAseq/**
      - Elbow plot - genes.png
      **Input/**
        - LFC_AllGenes - p-value.csv
        - input-values.csv
      - STRING-2Clusters.png
      - PCA-2Clusters.png
      - t-SNE-2Clusters.png
      - clustered_results.csv
    **Clusterización - Transcripción_fisiológica/**
      - Elbow.png
      - PCA.png
      **Input/**
        - input-values.csv
      - Scree.png
      - STRING.png
      - clustered_results.csv
    **Clusterización - Proteómica_cáncer/**
      - Elbow.png
      - PCA.png
      **Input/**
        - input-values.csv
      - t-SNE.png
      - Scree.png
      - STRING.png
      - clustered_results.csv
**Datos_Proteómicos/**
  - Análisis_STRING.py
  - Expresión_Proteínas.R
  **Datos/**
    - Expresión proteica_IHC - cáncer.csv
    - Expresión proteica_IHC - fisiológica.csv
**Datos_Genómicos/**
  - Gráfico_Alteraciones_Genéticas_ENPP2.py
  - Co-ocurrencia_Genómica.R
  - Gráfico_Análisis_Supervivencia_Genómica.py
  **Datos/**
    **Alteraciones Genómicas/**
      - Alteraciones Genómicas - Todos los genes.csv
      - Alteraciones Genómicas - ENPP2.xlsx
      - RawData - Tipos de alteraciones.xlsx
    **Co-ocurrencia/**
      - Co-ocurrencia_TodosLosGenes.csv
  - Gráfico_Alteraciones_Genéticas_Genes.py
