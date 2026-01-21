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

library(msigdbr)  #install.packages("msigdbr")
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)

## example
load('./sc_m1_from_macro.Rdata')

dat <- as.matrix(m1.from_macro@assays$RNA@counts)

geneset <- mouse.symbol.list

gsva_mat <- gsva(expr=dat, 
               gset.idx.list=geneset, 
               kcdf="Poisson" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
               verbose=T, 
               parallel.sz = parallel::detectCores())