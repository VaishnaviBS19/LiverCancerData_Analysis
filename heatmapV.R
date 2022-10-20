getwd()
# heatmap visualization of data, for this we need to download and load package  "complexheatmap"
#install.packages("ComplexHeatmap")
library(ComplexHeatmap)

# using 'source' we can use data from another r file.
data=source("D:\\Vaishnavi\\MSc sem 3\\Cancer_genomics\\CG_R\\dNormalization.r")

print(LogC)
library(matrixStats)

# calculating variance using logC
variance = apply(LogC,1,varDiff)
variance = sort(variance,decreasing = T)
Top50 = variance[1:50] # taking first 50 records to generate heatmap
zMat = ZScore[names(Top50),]
print(zMat)

# generate heatmap
Heatmap(zMat)






