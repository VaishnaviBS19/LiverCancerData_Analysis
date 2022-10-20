# SSGSEA  (Single Sample Geneset Enrichment Analysis)
library(matrixStats)
library(circlize)
library(ComplexHeatmap)
library(data.table)

CancerData=read.csv("GSE185853.csv",header = T,row.names = 1)

gene_sets = CancerData

ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

data = readRDS("logC.rds")
data[1:5,1:5]
data = as.matrix(data)
gene_set = read.csv("GSE185853.csv")
head(gene_set)
gene_sets = as.list(as.data.frame(gene_set))
print("genes set ready")
system.time(assign('res', ssgsea(data, gene_sets, scale = TRUE, norm = FALSE)))
res1 = t(res)
head(res1)

#zscore the ssgsea output for comparative analysis
mat = (res - rowMeans(res))/(rowSds(as.matrix(res)))[row(res)]

#read the info file for Heatmap annotations
info = read.csv("sampleInfoCGGA.csv")
rownames(info) = info[,1]
info = info[,-1]

#order both the objects for sample alignment. If the number of samples vary in your data, please subset the data frames and then order them.

mat = mat[, order(colnames(mat))]
info = info[order(rownames(info)), ]

#move on for further analysis, only if the statement is TRUE
identical(rownames(info), colnames(mat))                     

