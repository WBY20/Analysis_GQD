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
    library(cluster)
    library(stringr)
    library(GOSemSim)
    # library(DOSE)
    library(org.Hs.eg.db)
    library(R.matlab)
    library(igraph)
    library(poibin)
}
)

load('./macro1229.Rdata')
load('./mono1229.Rdata')

### generate file for SCENIC
write.table(t(as.matrix(dt@assays$RNA@data)), file = './TF/mono/expression.tsv', sep = '\t')
write.table(t(as.matrix(dtt@assays$RNA@data)), file = './TF/macro/expression.tsv', sep = '\t')

### Running in Linux CMD line

# read -p "Please enter user name:" user
# echo $user welcome

# for sample in TF
# do

#   # grn
#   pyscenic grn --num_workers 18 -o ./$sample/expr_mat.tsv ./$sample/expression.tsv ./allTFs_hg38.txt

#   # ctx
#   pyscenic ctx ./$sample/expr_mat.tsv ./hg19-tss-centered-5kb-7species.mc9nr.feather ./hg19-tss-centered-10kb-7species.mc$

#   # aucell
#   pyscenic aucell ./$sample/expression.tsv ./$sample/regulons.csv -o ./$sample/auc_mtx.csv --num_workers 24

#   echo $sample has finished

# done

# output = Work Finished

# echo $output

### Finding TFs
library(limma)
library(edgeR)

tf1 = read.table('./TF/mono/regulons.csv', sep=',')
tf2 = read.table('./TF/macro/regulons.csv', sep=',')

### example
m2.gqd = intersect(rownames(auc_mtx), colnames(m2.macro[,which(m2.macro$group == 'GQD')]))
m2.bbp = intersect(rownames(auc_mtx), colnames(m2.macro[,which(m2.macro$group == 'BBP')]))
m2.t2dm = intersect(rownames(auc_mtx), colnames(m2.macro[,which(m2.macro$group == 'D')]))

auc_mtx = read.csv('./TF/mono/auc_mtx.csv', row.names = 1)
auc_mtx = auc_mtx[union(m2.gqd, m2.t2dm), ]

auc.analysis <- function(auc_mtx){
auc_mtx_raw <- auc_mtx[,-which(colnames(auc_mtx) %in% c("ident"))]
auc_mtx_raw <- t(auc_mtx_raw)
d <- DGEList(auc_mtx_raw)
d <- calcNormFactors(d)
groups <- as.character(auc_mtx$ident)

all.idents <- unique(groups)
result.list <- list()
for(i in all.idents){
    groups.temp <- groups
    groups.temp[which(groups == i)] <- "Case"
    groups.temp[-which(groups == i)] <- "Norm"
    mm <- model.matrix(~0 + groups.temp)
    y <- voom(d, mm)
    fit <- lmFit(y, mm)
    contr <- makeContrasts(groups.tempCase-groups.tempNorm, levels = colnames(coef(fit)))
   tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    p.val.results <- apply(auc_mtx_raw,1,function(x){p.val <- t.test(x[which(groups == i)],x[-which(groups == i)],alternative = "greater"); p.val <- p.val$p.value;return(p.val)})
    top.table$t.test.pval <- p.val.results[rownames(top.table)]
    top.table$t.test.pval.adjust <- p.adjust(top.table$t.test.pval)
    result.list[[i]] <- top.table
}    
return(result.list)
}

auc.analysis.results <- auc.analysis(auc_mtx)

TF2 = auc.analysis.results[['GQD']][auc.analysis.results[['GQD']]$adj.P.Val < 0.05, ]
rownames(TF2) = gsub('\\...','',rownames(TF2))
