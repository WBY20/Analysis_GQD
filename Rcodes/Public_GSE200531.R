options(warn = -1)
suppressPackageStartupMessages(
{
    library(Seurat)
    library(cowplot)
    library(Matrix)
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
    library(CytoTRACE)
    library(tidyr)
    library(ggnewscale)
    library(patchwork)
    library(tidydr)
    library(ComplexHeatmap)
    library(ggrepel)
    library(factoextra)
    library(future)
    library(future.apply)
    library(fpc)
    library(clusterProfiler)
    library(cluster)
    library(stringr)
    library(limma)
    library(edgeR)
    library(GOSemSim)
    # library(DOSE)
    library(org.Hs.eg.db)
    library(R.matlab)
    library(igraph)
    library(poibin)
}
)



packageVersion('Seurat')

human2mouse <- read.csv('./Immune/human2mouse.csv')
gconvert <- function(input){
    ss <- bitr(input, fromType = "SYMBOL",toType = c( "ENSEMBL"),OrgDb = org.Hs.eg.db)
    ee <- human2mouse[which(human2mouse$Gene.stable.ID %in% ss$ENSEMBL),]$Mouse.gene.name
    gg <- unique(ee)
    return(gg)
}



T2DM.data = Read10X(data.dir = "./GSE200531/T2DM/")

T2DM <- CreateSeuratObject(counts = T2DM.data, project = "T2DM", min.cells = 3, min.features = 200)

T2DM[["percent.mt"]] <- PercentageFeatureSet(T2DM, pattern = "^mt-")

VlnPlot(T2DM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = -1)

T2DM <- subset(T2DM, subset = percent.mt < 25)

VlnPlot(T2DM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = -1)

T2DM <- subset(T2DM, subset = percent.mt < 10)

T2DM

T2DM = NormalizeData(T2DM, normalization.method = "LogNormalize", scale.factor = 10000)

T2DM <- FindVariableFeatures(T2DM, selection.method = "vst", nfeatures = 2000)

T2DM <-  ScaleData(T2DM,verbose = F) %>% RunPCA(npcs = 50,verbose = F) %>%
    RunUMAP(dims = 1:20,verbose = F) %>% FindNeighbors(dims = 1:20,verbose = F) %>% FindClusters(resolution = 0.2,verbose = F)
T2DM <- RunTSNE(T2DM, dims = 1:20)

T2DM$site = 'T2D'

DimPlot(T2DM, reduction = 'umap', label = T)

FeaturePlot(T2DM, features = gconvert(c('VIM', 'EPCAM', 'PTPRC')))



Im.markers = FindAllMarkers(T2DM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log2(1.5))

Im.markers[Im.markers$cluster == '6', ]

new.cluster.ids <- c("Non Immune cells", "Non Immune cells", "B cells", "Non Immune cells", "Non Immune cells", "Fibroblasts", 
                     "Non Immune cells", "Non Immune cells", "Non Immune cells", 'Macrophages', 'CD8+ T', 'CD4+ T')
names(new.cluster.ids) <- levels(T2DM)
T2DM <- RenameIdents(T2DM, new.cluster.ids)
DimPlot(T2DM, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

save(T2DM, file = './GSE200531.Rdata')

load('./GSE200531.Rdata')

options(repr.plot.height =6, repr.plot.width = 6)
(DimPlot(T2DM, reduction = 'umap', label = T, label.size = 6) + theme_dr() + NoLegend() 
 + theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank()) 
)

T2DM$GQD = colMeans(T2DM@assays$RNA@data[intersect(rownames(T2DM),gconvert(gqd.symbol$SYMBOL)),])

VlnPlot(T2DM, features = 'GQD', sort = T,  pt.size = -1) + NoLegend()

T2DM = AddModuleScore(T2DM, features = list(gconvert(gqd.symbol$SYMBOL)), name = 'GQD', assay = 'RNA')

options(repr.plot.height =6, repr.plot.width = 6)
(VlnPlot(T2DM, features = 'GQD1', sort = T,  pt.size = -1, idents = c('Macrophages', 'CD4+ T', 'CD8+ T', 'B cells')) + theme_bw() + labs(title = '', x = '')
 + theme(panel.grid = element_blank()) + theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20, angle = 45, hjust = 1))
)+ NoLegend() 

options(repr.plot.height =4, repr.plot.width = 6)
(VlnPlot(T2DM, features = 'GQD1', sort = T,  pt.size = -1, idents = c('Macrophages', 'CD4+ T', 'CD8+ T', 'B cells')) + theme_bw() + labs(title = '', x = '')
 + theme(panel.grid = element_blank()) + theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20, angle = 45, hjust = 1))
)+ NoLegend() 

options(repr.plot.height =4, repr.plot.width = 8)
(VlnPlot(T2DM, features = 'GQD1', sort = T,  pt.size = -1, idents = c('Macrophages', 'CD4+ T', 'CD8+ T', 'B cells')) + theme_bw() + labs(title = '', x = '')
 + theme(panel.grid = element_blank()) + theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20, angle = 30, hjust = 1))
)+ NoLegend() 

