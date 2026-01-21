options(warn = -1)
suppressPackageStartupMessages(
{
    # library(Seurat)
    library(cowplot)
    library(Matrix)
    library(dplyr)
    # library(monocle)
    library(reshape2)
    library(pheatmap)
    library(ggbiplot)
    library(VGAM)
    library(reshape2)
    library(ggplot2)
    library(mgcv)
    library(latex2exp)
    library(ggpubr)
    # library(CytoTRACE)
    library(tidyr)
    # library(ggnewscale)
    library(patchwork)
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
    # library(GOSemSim)
    # library(DOSE)
    library(org.Hs.eg.db)
    library(R.matlab)
    library(igraph)
    library(poibin)
}
)

library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont=c("BP"))

library(parallel)
parallel.cores <- detectCores()  
# cl <- makeCluster(parallel.cores-2)
cl <- makeCluster(36)

clusterEvalQ(cl, { library('GOSemSim') })  =
clusterExport(cl, c('cmd.target', 'ref.gene', 'hsGO'), envir = environment())  

gqd.sim = parLapply(cl, cmd.target, function(x){
    gene = x
    score = clusterSim(gene, ref.gene$ENTREZID, semData=hsGO, measure="Wang", combine="BMA")
    return(score)
})

stopCluster(cl)
