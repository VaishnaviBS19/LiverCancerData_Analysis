#single cell 
#setwd("D:/Vaishnavi/MSc sem 3/Cancer_genomics/CG_R")
#install.packages('Seurat')
#install.packages('dplyr')
#install.packages('patchwork')

library(Seurat)
library(dplyr)
library(patchwork)

pdata = Read10X(data.dir = "singleCellAnalysis")  # 10X means - 10 times each sequence from your data is read
print(pdata)


# creating object of seurat
ObData<- CreateSeuratObject(counts = pdata, min.cells = 3, min.features = 200)
print(ObData)

pdata[1:50, 1:10]

ObData[["percent.mt"]] = PercentageFeatureSet(ObData, pattern = "^MT-")
head(ObData@meta.data)


# Visualize QC metrics as a violin plot
VlnPlot(ObData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 = FeatureScatter(ObData, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(ObData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

ObData = subset(ObData, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
ObData

ObData = NormalizeData(ObData)
ObData = FindVariableFeatures(ObData, selection.method = "vst", nfeatures = 2000)

top10 = head(VariableFeatures(ObData), 10)
top10

# plot variable features with and without labels
plot1 = VariableFeaturePlot(ObData)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1


all.genes = rownames(ObData)
ObData = ScaleData(ObData, features = all.genes)

ObData@assays$RNA@scale.data[1:50, 1:5]

ObData = FindNeighbors(ObData, dims = 1:10)

ObData = FindClusters(ObData, resolution = 0.5)

head(ObData@meta.data)

ObData = RunUMAP(ObData, dims = 1:10)

DimPlot(ObData, reduction = "umap")


pbmc.markers = FindAllMarkers(ObData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

a = pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
a

genes = a %>% pull(gene)
genes

FeaturePlot(pbmc, features = genes[1:2])

FeaturePlot(pbmc, features = genes[1:2], cols = c("white", "red"))

