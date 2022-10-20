# DEG analysis(Analysis of differentially expressed genes) and plotting volcano plot

data = source("heatmapV.r")
# converting data of logC values into matrix
Lmatrix = matrix(NA,ncol = 4,nrow = nrow(LogC))
print(Lmatrix)
rownames(Lmatrix)=rownames(LogC)
colnames(Lmatrix)=c('LiverTumor','Control','pvalue','log2FC')    #log2FC = log2foldchange

for (i in 1:nrow(LogC)) {
  vector1 = as.numeric(LogC[i,1:3])
  vector2 = as.numeric(LogC[i,4:6])
  
  test1 = t.test(vector1,vector2,paired = F,alternative = "two.sided")
  Lmatrix[i,1]=test1$estimate[[1]]
  Lmatrix[i,2]=test1$estimate[[2]]
  Lmatrix[i,3]=test1$p.value
  Lmatrix[i,4]=Lmatrix[i,1]-Lmatrix[i,2]
  
}

Lmatrix = as.data.frame(Lmatrix)
num=which(is.nan(Lmatrix$pvalue))
Lmatrix[num,'pvalue'] = 1

# save Lmatrix to rds format
saveRDS(Lmatrix,file = "Lmatrix.rds")


# Visualizing of differentially expressed genes using Volcano plot  using package 'EnhancedVolcano'

library(EnhancedVolcano)
vplot = EnhancedVolcano(Lmatrix,lab = rownames(Lmatrix),x = 'log2FC', y = 'pvalue' )
print(vplot)
