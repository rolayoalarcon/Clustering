# Clustering
## Scripts:
Hemos utilizado 3 scripts principales. Primeramente existe **get_bits_heavy.py**: Con este script realizamos el BLASTP para obtener los bitscores descritos en el trabajo. Tambien se arma la matriz simetrica _pesada_ de bitscores, utilizando el puntaje mas alto obtenido para cada par de proteinas. Por ultimo se llama al script **Cluster_Protein_Analysis.R** para el analisis.

En **Cluster_Protein_Analysis.R** realizamos en analisis de clustering para la matriz pesada (ver reporte escrito). En este script se hace el clustering jerarquico con los distintos metodos disponibles. Se crean los dendogramas, heatmaps y scatterplots correspondientes.
En **Cluster_Protein_Analysis_Light.R** se lleva a cabo el mismo analisis para la matriz ligera

## Resultados
En **NormMatrix** se encuentran las figuras generadas por **Cluster_Protein_Analysis.R** para la matriz pesada, normalizada con el valor maximo de la matriz
En **NormRow** se encuentran las figuras generadas por **Cluster_Protein_Analysis.R** para la matriz pesada, normalizada con el valor maximo de cada renglon

En **LightNormMatrix** se encuentran las figuras generadas por **Cluster_Protein_Analysis_Light.R** para la matriz ligera, normalizada con el valor maximo de la matriz
En **NormRow** se encuentran las figuras generadas por **Cluster_Protein_Analysis_Light.R** para la matriz ligera, normalizada con el valor maximo de cada renglon
