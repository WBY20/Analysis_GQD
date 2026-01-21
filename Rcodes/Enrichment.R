library(clusterProfiler)

library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(org.Rn.eg.db)

## example
load('./DEGs.Rdata') ## load DEGs

diff.gene = bitr(diff.gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

### KEGG
kegg = enrichKEGG(gene = diff.gene$ENTREZID, keyType = "kegg", organism = 'hsa', pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff  = 0.05)
kegg = kegg@result

kegg = kegg[kegg$p.adjust < 0.05, ]

### GO
go = enrichGO(OrgDb="org.Hs.eg.db", diff.gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 0.05, readable= TRUE)
go = go@result