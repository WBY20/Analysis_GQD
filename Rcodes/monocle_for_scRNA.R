options(warn = -1)
suppressPackageStartupMessages(
{
    library(Seurat)
    library(cowplot)
    library(Matrix)
    library(clusterProfiler)
    library(dplyr)
    library(monocle)
    library(reshape2)
    library(ggnewscale)
    library(pheatmap)
    library(ggbiplot)
    library(VGAM)
    library(reshape2)
    library(ggplot2)
    library(mgcv)
    library(latex2exp)
    library(ggpubr)
    library(tidyr)
    library(ggnewscale)
    library(patchwork)
    library(ComplexHeatmap)
    library(ggrepel)
    library(factoextra)
    library(cluster)
    library(stringr)
    library(GOSemSim)
    library(DOSE)
    library(org.Hs.eg.db)
    library(R.matlab)
    library(igraph)
}
)

load('./macro1229.Rdata')

options(repr.plot.height=6, repr.plot.width=6)
DimPlot(dtt, reduction = 'umap', label = T, repel = T)

dtt$final_celltype = 'M1 macro'
dtt$final_celltype[Idents(dtt) == 'M2 macro'] = 'M2 macro'
dtt$final_celltype[Idents(dtt) == 'M0 macro'] = 'M0 macro'

dtt$macrotype = Idents(dtt)

m2 = subset(dtt, idents = 'M2 macro')



m2

m2 = FindVariableFeatures(m2, selection.method = "vst", nfeatures = 2000)
# data <- ScaleData(data, features = all.genes)
m2 <- ScaleData(m2)

m2 <- RunPCA(m2, verbose = FALSE, features = VariableFeatures(object = m2), npcs = 50)
m2 <- FindNeighbors(m2, dims = 1:20, verbose = FALSE)
m2 <- FindClusters(m2, verbose = FALSE, resolution = 1)
m2 = RunUMAP(m2, dims = 1:20)

options(repr.plot.height=12, repr.plot.width=12)
FeaturePlot(m2, features = c('M1.score1', 'M2.score1', 'Mrc1', 'M0', 'M1', 'M2', 'nFeature_RNA', 'Tnf', 'Arg1'), 
            min.cutoff = 0, label = T, ncol = 3)

options(repr.plot.height=12, repr.plot.width=12)
VlnPlot(m2, features = c('M1.score1', 'M2.score1', 'Mrc1', 'M0', 'M1', 'M2', 'nFeature_RNA', 'Tnf', 'Arg1'), sort = T, ncol = 3)



library(CytoTRACE)

example.data = data.frame(as.matrix(m2@assays$RNA@counts))

results.test1 <- CytoTRACE(example.data, ncores = 8, subsamplesize = 1000)

m2$stemness = results.test1$CytoTRACE

options(repr.plot.height=6, repr.plot.width=6)
FeaturePlot(m2, features = 'stemness')
VlnPlot(m2, features = 'stemness', sort = T)



count_matrix <- m2@assays$RNA@counts
pd <- m2@meta.data
fd <- m2@assays$RNA@meta.features
pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)
options(warn = -1)
# mono.dt.epi <- newCellDataSet(as.matrix(count_matrix),
#                 phenoData = pd,
#                 featureData = fd,
#                 lowerDetectionLimit = 0.1,
#                 expressionFamily=tobit(Lower = 0.1))
# rpc_matrix <- relative2abs(mono.dt.epi, method = "num_genes")

mono.dt.epi <- newCellDataSet(as(as.matrix(count_matrix), "sparseMatrix"),
                phenoData = pd,
                featureData = fd,
                lowerDetectionLimit = 0.5,
                expressionFamily = negbinomial.size())
fData(mono.dt.epi)$gene_short_name <- rownames(fData(mono.dt.epi))
mono.dt.epi <- estimateSizeFactors(mono.dt.epi)
mono.dt.epi <- estimateDispersions(mono.dt.epi)
mono.dt.epi <- detectGenes(mono.dt.epi, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(mono.dt.epi),num_cells_expressed >= 10)) ### expressed in 10 cells at least 

# # valid_cells <- row.names(subset(pData(mono.dt.epi),
# #             num_genes_expressed >= quantile(dtt$nFeature_RNA)['50%'] &
# #             nCount_RNA >= quantile(dtt$nCount_RNA)['50%']))
# valid_cells <- row.names(subset(pData(mono.dt.epi),
#             num_genes_expressed > 300 &
#             nCount_RNA > 5000))
# length(valid_cells)
# # length(data.temp.1$sample)
# mono.dt.epi <- mono.dt.epi[,valid_cells]

pData(mono.dt.epi)$Total_mRNAs <- Matrix::colSums(exprs(mono.dt.epi))
mono.dt.epi <- mono.dt.epi[,pData(mono.dt.epi)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(mono.dt.epi)$Total_mRNAs)) +
            2*sd(log10(pData(mono.dt.epi)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(mono.dt.epi)$Total_mRNAs)) -
            2*sd(log10(pData(mono.dt.epi)$Total_mRNAs)))
qplot(Total_mRNAs, data = pData(mono.dt.epi), geom ="density") +
geom_vline(xintercept = lower_bound) +
geom_vline(xintercept = upper_bound)

mono.dt.epi <- mono.dt.epi[,pData(mono.dt.epi)$Total_mRNAs > lower_bound &
      pData(mono.dt.epi)$Total_mRNAs < upper_bound]
mono.dt.epi <- detectGenes(mono.dt.epi, min_expr = 0.1)

dim(mono.dt.epi)

### variable genes as 500 variable genes in each cell cluster
options(repr.plot.width = 8, repr.plot.height = 4)
dtt <-  FindVariableFeatures(dtt, nfeatures = 1000, assay = 'RNA')
order.genes <- dtt@assays$RNA@var.features
length(order.genes)
mono.dt.epi <- setOrderingFilter(mono.dt.epi, order.genes)
plot_ordering_genes(mono.dt.epi)
mono.dt.epi <- reduceDimension(mono.dt.epi, max_components = 2, method = 'DDRTree')
mono.dt.epi <- orderCells(mono.dt.epi)
plot_cell_trajectory(mono.dt.epi, color_by = "seurat_clusters")
# plot_cell_trajectory(mono.dt.epi, color_by = "group")
# plot_cell_trajectory(mono.dt.epi, color_by = "macrotype")

plot_cell_trajectory(mono.dt.epi, color_by = "M1")
plot_cell_trajectory(mono.dt.epi, color_by = "stemness")
plot_cell_trajectory(mono.dt.epi, color_by = "State")





plot_cell_trajectory(mono.dt.epi, color_by = "group")

mono.dt.epi <- orderCells(mono.dt.epi, root_state = 1)
plot_cell_trajectory(mono.dt.epi, color_by = "Pseudotime")

p1 = plot_cell_trajectory(mono.dt.epi, color_by = "M2")
p2 = plot_cell_trajectory(mono.dt.epi, color_by = "Pseudotime")
options(repr.plot.width = 8, repr.plot.height = 4)
plot_grid(plotlist = list(p1, p2), ncol = 2)

colnames(m2@meta.data)







pseudo <- data.frame(row.names = rownames(pData(mono.dt.epi)))
pseudo$pseudotime <- pData(mono.dt.epi)$Pseudotime
pseudo$m21 <- m2@meta.data[rownames(pData(mono.dt.epi)), ]$M2
pseudo$m22 <- m2@meta.data[rownames(pData(mono.dt.epi)), ]$M2.score1
pseudo$m11 <- m2@meta.data[rownames(pData(mono.dt.epi)), ]$M1
pseudo$m12 <- m2@meta.data[rownames(pData(mono.dt.epi)), ]$M1.score1
pseudo$group = m2@meta.data[rownames(pData(mono.dt.epi)), ]$group



pp = pseudo

pp$pseudotime = (pp$pseudotime - min(pp$pseudotime)) / (max(pp$pseudotime) - min(pp$pseudotime))

options(repr.plot.width = 8, repr.plot.height = 4)
(ggplot(pp, aes(x= pseudotime, y= m21, color = group)) 
  +labs(title = 'M2 score in different Groups',x= 'Normalized progression pseudotime', y = 'M2 score')+theme(plot.title = element_text(hjust = 0.5,size = 25))
  + geom_smooth(level = 0.95, method = 'loess',show.legend= FALSE, se = F) 
  + geom_smooth(fill = NA, method = 'loess')
# +scale_color_manual(name = 'Pathway',values = path.color,breaks = path.list, labels=path.label)
  +theme(axis.text = element_text(size=25),axis.title = element_text(size=25))
  +theme(legend.text = element_text(size=25),legend.title =element_text(size=25),legend.key = element_rect(fill='white'))
  +theme(panel.background = element_blank(),axis.line = element_line()))







sub.mono.dt.epi = rownames(pData(mono.dt.epi)[pData(mono.dt.epi)$State %in% c(1,3), ])
pseudo <- data.frame(row.names = sub.mono.dt.epi)
pseudo$pseudotime <- pData(mono.dt.epi)[pData(mono.dt.epi)$State %in% c(1,3), ]$Pseudotime
pseudo$m21 <- m2@meta.data[sub.mono.dt.epi, ]$M2
pseudo$m22 <- m2@meta.data[sub.mono.dt.epi, ]$M2.score1
pseudo$m11 <- m2@meta.data[sub.mono.dt.epi, ]$M1
pseudo$m12 <- m2@meta.data[sub.mono.dt.epi, ]$M1.score1
pseudo$group = m2@meta.data[sub.mono.dt.epi, ]$group



pp = pseudo

pp$pseudotime = (pp$pseudotime - min(pp$pseudotime)) / (max(pp$pseudotime) - min(pp$pseudotime))

options(repr.plot.width = 8, repr.plot.height = 4)
(ggplot(pp, aes(x= pseudotime, y= m11, color = group)) 
  +labs(title = 'M2 score in different Groups',x= 'Normalized progression pseudotime', y = 'M2 score')+theme(plot.title = element_text(hjust = 0.5,size = 35))
  + geom_smooth(level = 0.95, method = 'loess',show.legend= FALSE, se = F) 
  + geom_smooth(fill = NA, method = 'loess')
# +scale_color_manual(name = 'Pathway',values = path.color,breaks = path.list, labels=path.label)
  +theme(axis.text = element_text(size=25),axis.title = element_text(size=25))
  +theme(legend.text = element_text(size=25),legend.title =element_text(size=25),legend.key = element_rect(fill='white'))
  +theme(panel.background = element_blank(),axis.line = element_line()))


