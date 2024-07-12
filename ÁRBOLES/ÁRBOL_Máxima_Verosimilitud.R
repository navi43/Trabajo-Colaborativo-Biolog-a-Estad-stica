# ÁRBOL CON MÁXIMA VEROSIMILITUD------------------------------------------------

# Cargar las librerías necesarias
library(ape)
library(phangorn)

# Cargar el archivo FASTA y el archivo CSV con anotaciones
file_path <- "secuencia_alineada.fasta"  # Ruta al archivo FASTA con secuencias de pasiflora
csv_path <- "secuencia_alineada.csv"     # Ruta al archivo CSV con anotaciones

dna <- read.dna(file_path, format = "fasta")
annot <- read.csv(csv_path, header = TRUE, row.names = 1)

# Convertir las secuencias a formato phyDat para análisis filogenético
dna_phy <- as.phyDat(dna)

# Crear un árbol inicial usando Neighbor Joining
tre_ini <- nj(dist.dna(dna, model = "TN93"))
# Se crea un árbol filogenético inicial utilizando el método de Neighbor Joining
# basado en las distancias de ADN calculadas con el modelo TN93.

# Estimar la verosimilitud del árbol inicial
fit_ini <- pml(tre_ini, dna_phy, k = 4)
# Se estima la verosimilitud utilizando el método de Máxima Verosimilitud 
# con 4 categorías de gamma.

# Optimizar el árbol inicial
fit <- optim.pml(fit_ini, optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE)
# tasas de sustitución de nucleótidos (optBf)
# las frecuencias de bases (optQ)
# la forma de la distribución gamma (optGamma)


# Mostrar los resultados de la optimización
print(fit)

# Comparar los modelos inicial y optimizado usando AIC (criterio de información de Akaike)
aic_ini <- AIC(fit_ini)
aic_opt <- AIC(fit)
cat("AIC inicial:", aic_ini, "\nAIC optimizado:", aic_opt, "\n")
# Un AIC menor indica un mejor ajuste del modelo.

# Enraizar y organizar el árbol optimizado
tre_opt <- root(fit$tree, 1)
tre_opt <- ladderize(tre_opt)

# Definir una paleta de colores pasteles
colores <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# Crear un gráfico del árbol optimizado
windows()
par(mar = c(5, 4, 4, 2) + 0.1) # Ajustar márgenes
plot(tre_opt, show.tip = FALSE, edge.width = 2, cex = 0.8) # Mostrar etiquetas de las especies
#show.tip = TRUE -> Para poder observar a que tipo de especie de passiflora nos referimos

# Obtener nombres únicos de las especies
unique_species <- unique(annot$ID)

# Asignar colores a cada especie
species_colors <- colores[match(annot$ID, unique_species)]

# Añadir etiquetas de las especies con colores diferentes
tiplabels(annot$ID, bg = species_colors, col = "black", cex = 0.7)

# Añadir el eje temporal
axisPhylo()

# Ajustar título del gráfico
title("Árbol de Máxima Verosimilitud de Pasiflora", cex.main = 1.2)
