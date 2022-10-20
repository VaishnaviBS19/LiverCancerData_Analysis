setwd("D:\\Vaishnavi\\MSc sem 3\\Cancer_genomics\\CG_R")
getwd()

# read liver cancer input file
CancerData=read.csv("GSE185853.csv",header = T,row.names = 1)
print(CancerData)
summary(CancerData)

#calculate count per matrix(CPM)
CPM = CancerData
for(i in 1:ncol(CancerData)){
  CPM[,i]=(CancerData[,i]/sum(CancerData[,i]))*1000000
}

print(CPM)

#log of CPM
LogC=log2(CPM+1)
print(LogC)
#save it in RDS form
saveRDS(LogC,file = "logC.rds")

# z score calculation using log of cpm
# for calculating z score download and load package matrixStats
library(matrixStats)
ZScore = (LogC - rowMeans(LogC))/rowSds(as.matrix(LogC))[row(LogC)]
print(ZScore)

# Now our data is NORMALIZED




