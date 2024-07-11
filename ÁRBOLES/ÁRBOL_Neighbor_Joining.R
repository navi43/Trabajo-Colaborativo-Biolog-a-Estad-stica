# RECORDAR: Usamos las secuencias que elegimos, pero para que funcione el código
# primero tenemos que alinearlas en un programa llamado MAFFT.

# ÁRBOL CON MÉTODO Neighboor Join-----------------------------------------------
library(adegenet)
file_path <- "secuencia_alineada.fasta"
dna <- fasta2DNAbin(file=file_path) 
# convierte el archivo FASTA en un objeto de clase DNAbin que puede ser usado 
# para análisis filogenéticos.
dna
csv_path<- "secuencia_alineada.csv"
annot <- read.csv(csv_path, header=TRUE,row.names = 1)
# lee un archivo csv que contiene anotaciones relacionadas con las secuencias de
# ADN
annot

library(ape)
D <- dist.dna(dna, model = "TN93")

length(D)
temp <- as.data.frame(as.matrix(D)) # Convierte la matriz de distancias en un 
                                    # data frame para su visualización.

print(D)

#Visualización de la matriz de distancia 
print(temp)

# Cargar el paquete ade4
library(ade4)
windows()
table.paint(temp, cleg=0, clabel.row=0.5, clabel.col=0.5)
# la matriz de distancias se visualiza como una tabla pintada, donde los valores
# de la matriz se representan con colores.

#Construyendo árbol utilizando NJ:
library(ape)
tree <- nj(D) # Método de NJ para construir un árbol filogenético basado en la 
              # matriz de distancias.
class(tree) 
tree <- ladderize(tree) # se reorganiza el árbol
tree

plot(tree, cex = 0.6) # cex: escala de caracteres
title("Un simple Arbol NJ")

#Arbol no enraizado
library(ape)
library(adegenet)
windows(width = 8, height = 5)
plot(tree, show.tip = FALSE) # sin mostrar el nombre de las secuencias
title("Arbol NJ no enraizado")
myPal <- colorRampPalette(c("red", "yellow", "green", "blue"))
tiplabels(annot$ID, bg = myPal(factor(annot$ID)), cex = 0.5)
legend("topright", fill = myPal(length(unique(annot$ID))), 
       legend = unique(annot$ID), ncol = 2)

#Árbol enraizado con una de las anotaciones
library(ape)
library(adegenet)
windows(width = 8, height = 5)
tree2 <- root(tree, out = 1)
tree2 <- ladderize(tree2)
windows(width = 8, height = 5)
plot(tree2, show.tip=FALSE, edge.width=1)
title("Arbol enraizado NJ")
myPal <- colorRampPalette(c("red", "yellow", "green", "blue"))
tiplabels(annot$ID, bg = myPal(factor(annot$ID)), cex = 0.5)
legend("bottomright", fill = myPal(length(unique(annot$ID))), 
       legend = unique(annot$ID), ncol = 2)
axisPhylo() # Añade una leyenda y un eje filogenético (axisPhylo)