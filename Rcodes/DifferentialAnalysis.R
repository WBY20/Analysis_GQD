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
    library(ComplexHeatmap)
    library(ggrepel)
    library(factoextra)
    library(future)
    library(future.apply)
    library(fpc)
    library(clusterProfiler)
    library(KEGGREST)
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

## example : AA
dt = read.xlsx('./2024-054-C-01_202404180515_bioNormConc.xlsx')

normal = colnames(dt)[grep('C', colnames(dt))]
t2dm = colnames(dt)[grep('D', colnames(dt))]
gqd = colnames(dt)[grep('GM', colnames(dt))]
bbp = colnames(dt)[grep('BBP', colnames(dt))]

comparison = list('D vs C'=c('D','C'), 'GQD vs D'=c('GM','D'), 'BBP vs D'=c('BBP','D'))
mapping = list('C'=normal, 'D'=t2dm, 'GM'=gqd, 'BBP'=bbp)

result = data.frame(meta=0, p=0, logfc=0, group=0, sp=0, s=0)
result = result[-1, ]
for(i in 1:length(comparison)){
    com = comparison[[i]]
    group1 = com[1]
    group2 = com[2]
    group1 = mapping[[group1]]
    group2 = mapping[[group2]]
    for(j in 1:dim(dt)[1]){
        meta = dt[j, 1]
        g1 = dt[j, group1]
        g2 = dt[j, group2]
        test = wilcox.test(as.numeric(g1), as.numeric(g2), paired = F) ### two-side
        logfc = log2(mean(as.numeric(g1))/mean(as.numeric(g2)))
        side = ifelse(logfc > 0, 'g', 'l')
        test2 = wilcox.test(as.numeric(g1), as.numeric(g2), paired = F, alternative = side) ### single-side
        result = rbind(result, data.frame(meta=meta, p=test$p.value, logfc=logfc, 
                                          group=names(comparison)[i], sp=test2$p.value, s=side))
    }
}

### adjustment
result$padj = p.adjust(result$p, method = "BH")
result$spadj = p.adjust(result$sp, method = "BH")

write.xlsx(result, './AA.xlsx')