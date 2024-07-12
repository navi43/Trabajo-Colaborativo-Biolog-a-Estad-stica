#ARBOL FILOGENÉTICO CON METODO UPGMA--------------------------------------------

# Instalar los paquetes si no están ya instalados
install.packages("ape")
install.packages("phangorn")

# Cargar los paquetes
library(ape)
library(phangorn)

# Leer el archivo de alineamiento en formato FASTA
alineamiento <- read.dna("secuencia_alineada.fasta", format = "fasta")

# Calcular la matriz de distancias usando el modelo K80
distancias <- dist.dna(alineamiento, model = "K80")
#Se calcula una matriz de distancias genéticas entre las secuencias de ADN
#Especifica el uso del modelo evolutivo Kimura 80 para calcular las distancias(model)

# Construir el árbol filogenético usando el método UPGMA
arbol_upgma <- upgma(distancias)

# Dibujar el árbol filogenético
plot(arbol_upgma, main = "Árbol Filogenético UPGMA")

