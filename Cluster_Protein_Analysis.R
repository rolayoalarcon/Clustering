library(methods)
library(cluster)
library(factoextra)
library(ctc)

#Importamos datos
BitMatrix <- read.table("BitMatrix.csv", sep = ',', header = T, row.names = 1)
clustermethods <- c("single", "average", "complete", "ward.D")

#Definir normalizaciones
NormByMatrix <- function(x){x <- max(BitMatrix)/x; return(x)}
NormByRow <- function(x){x <- max(x)/x; return(x)}

#Crear matrices normalizadas
NorMat <- apply(BitMatrix, 1, NormByMatrix)
NorRow <- apply(BitMatrix, 1, NormByRow)

#Exploración de datos
png("NormMatrix/Heatmap_ByMatrix.png",width = 777, height = 463)
heatmap(NorMat, symm = T)
dev.off()

png("NormRow/Heatmap_ByRow.png",width = 777, height = 463)
heatmap(NorRow, symm = T)
dev.off()


#Ahora si, los clusters 
AgCoMat <- c()
AgCoRow <- c()
for (met in clustermethods){
  
  clusMat <- hclust(as.dist(NorMat), method = met)
  clusRow <- hclust(as.dist(NorRow), method = met)
  
  #All for Matrix
    #Dendograma
  png(paste0(c("NormMatrix/"), met, c("_ClusterDendogram.png")), width = 777, height = 463)
  plot(clusMat, hang = -1)
  rect.hclust(clusMat, k=7, border = 3:4)
  dev.off()
  
    #Guardar newick
  newick_name <- file(paste0(c("NormMatrix/"), met, c("_Newick.txt")))
  writeLines(hc2Newick(clusMat),newick_name)
  close(newick_name)

    #Calcular coeficiente aglomerativo
  AgCoMat <- c(AgCoMat, met=coef.hclust(clusMat))
  
      #Scatter plot
  groupsMat <- cutree(clusMat, k=7)
  
  png(paste0(c("NormMatrix/"), met, c("_ScatterPlot.png")),width = 777, height = 463)
  p <- fviz_cluster(list(data = NorMat, cluster = groupsMat))
  print(p)
  dev.off()

  
  #All for Rows
  png(paste0(c("NormRow/"), met, c("_ClusterDendogram.png")), width = 777, height = 463)
  plot(clusRow, hang = -1)
  rect.hclust(clusRow, k=7, border = 3:4)
  dev.off()
  
  newick_name <- file(paste0(c("NormRow/"), met, c("_Newick.txt")))
  writeLines(hc2Newick(clusRow),newick_name)
  close(newick_name)
  
  #Calcular coeficiente aglomerativo
  AgCoRow <- c(AgCoRow, met=coef.hclust(clusRow))
  
  #Scatter plot
  groupsRow <- cutree(clusRow, k=7)
  
  png(paste0(c("NormRow/"), met, c("_ScatterPlot.png")), width = 777, height = 463)
  p <- fviz_cluster(list(data = NorRow, cluster = groupsRow))
  print(p)
  dev.off()
}
names(AgCoRow) <- clustermethods
names(AgCoMat) <- clustermethods
CompareACs <- rbind(AgCoMat, AgCoRow)
write.table(CompareACs, "AC_comparison.txt", sep = '\t')

