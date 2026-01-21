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

data <- readRDS("/data4/wby20/T2DM/5a.merge/dte.Mono.2.rds")

data@meta.data

m1 = read.csv('./M1.csv')
m2 = read.csv('./M2.csv')

dim(m1)
dim(m2)

options(repr.plot.height=6, repr.plot.width=6)
DimPlot(data, reduction = 'umap', label = T)

data



mono.cell = setdiff(levels(data), 'Prolif')

mono = subset(data, idents = mono.cell)

mono = FindVariableFeatures(mono, selection.method = "vst", nfeatures = 2000)
# data <- ScaleData(data, features = all.genes)
mono <- ScaleData(mono)

mono <- RunPCA(mono, verbose = FALSE, features = VariableFeatures(object = data), npcs = 50)
mono <- FindNeighbors(mono, dims = 1:20, verbose = FALSE)
mono <- FindClusters(mono, verbose = FALSE, resolution = 1)
mono = RunUMAP(mono, dims = 1:20)

options(repr.plot.height=6, repr.plot.width=8)
DimPlot(mono, reduction = 'umap', label = T)

table(mono$group)



# n = 20
# data = AddModuleScore(data, features = list(m1$Gene[1:20]), name = 'M1_top20')
# data = AddModuleScore(data, features = list(m2$Gene[1:20]), name = 'M2_top20')
# data = AddModuleScore(data, features = list(m1$Gene[1:10]), name = 'M1_top10')
# data = AddModuleScore(data, features = list(m2$Gene[1:10]), name = 'M2_top10')
# data = AddModuleScore(data, features = list(m1$Gene[1:5]), name = 'M1_top5')
# data = AddModuleScore(data, features = list(m2$Gene[1:5]), name = 'M2_top5')
mono = AddModuleScore(mono, features = list(m1$Gene), name = 'M1.score')
mono = AddModuleScore(mono, features = list(m2$Gene), name = 'M2.score')

options(repr.plot.height=8, repr.plot.width=8)
FeaturePlot(mono, features = c('M1.score1', 'M2.score1', 'M1', 'M2'),label = T, min.cutoff = 0)

options(repr.plot.height=4, repr.plot.width=8)
VlnPlot(mono, features = c('M1.score1', 'M2.score1'),pt.size = -1, sort = T)

options(repr.plot.height=8, repr.plot.width=12)
# FeaturePlot(data, features = c('Mrc1', 'Cd163', 'Mgl2', 'Cd86', 'Fcgr3', 'Cd80', 'Arg1', 'Nos2', 'Itgax'))
FeaturePlot(mono, features = c('Mrc1', 'Nkg7', 'Fabp4','S100a8','Cd209a', 'Xcr1'), label=T, ncol = 3)

# M0 and M1
options(repr.plot.height=12, repr.plot.width=12)
FeaturePlot(mono, features = c('Adgre1', 'Itgam', 'Cd68', 'Marco', 'Cd86', 'Cd80','Gpr18', 'Fpr2', 'Nos2'))
# FeaturePlot(data, features = c('M1.score1'), label=T, split.by = 'group')

# M2
options(repr.plot.height=12, repr.plot.width=12)
FeaturePlot(mono, features = c('Mrc1', 'Cd163', 'Itgam', 'Cd86', 'Tgm2', 'Egr2', 'Arg1', 'Myc', 'Mertk'))



count = as.data.frame(table(mono$group, Idents(mono)))

# m2.group = as.character(c(3, 10, 17, 20)) ## for 2000 variables
m2.group = as.character(c(3, 9, 12))
v1 = sum(count[which(count$Var2 %in% m2.group & count$Var1 == 'C'), ]$Freq)
v2 = sum(count[which(count$Var2 %in% m2.group & count$Var1 == 'D'), ]$Freq)
v3 = sum(count[which(count$Var2 %in% m2.group & count$Var1 == 'BBP'), ]$Freq)
v4 = sum(count[which(count$Var2 %in% m2.group & count$Var1 == 'GQD'), ]$Freq)
v1
v2
v3
v4

v1/3240
v2/5013
v3/4323
v4/3304

# m1.group = as.character(c(5, 18))
m1.group = as.character(c(5,13))
v1 = sum(count[which(count$Var2 %in% m1.group & count$Var1 == 'C'), ]$Freq)
v2 = sum(count[which(count$Var2 %in% m1.group & count$Var1 == 'D'), ]$Freq)
v3 = sum(count[which(count$Var2 %in% m1.group & count$Var1 == 'BBP'), ]$Freq)
v4 = sum(count[which(count$Var2 %in% m1.group & count$Var1 == 'GQD'), ]$Freq)
v1
v2
v3
v4

v1/3240
v2/5013
v3/4323
v4/3304

table(data$group)



options(repr.plot.height=12, repr.plot.width=12)
# FeaturePlot(data, features = c('Mrc1', 'Cd163', 'Mgl2', 'Cd86', 'Fcgr3', 'Cd80', 'Arg1', 'Nos2', 'Itgax'))
FeaturePlot(mono, features = c('Mrc1', 'Mertk', 'Myc', 'Fabp4','S100a8','Cd209a', 'Xcr1', 'Adgre1',
                               'Slc8a1', 'Cd163'), label=T, ncol = 3, min.cutoff = 0)

options(repr.plot.height=6, repr.plot.width=8)
DimPlot(mono, reduction = 'umap', label = T)



m1.group = c(5, 13)
m2a.group = c(3, 12)
m2c.group = c(9)
fabp4.group = c(4)
# s100a8.group = c(14)
cd209a.group = c(6, 10, 18)
m0.group = c(0, 1, 4, 7, 15)
# xcr1.group = c(11,14)
# cd24a.group = c(9)
# slc8a1.group = c(11)
# ccl5.group = c(4)
# ccl22.group = c(17)

dt = mono

new.cluster.ids = rep('Other Mono', length(levels(dt)))
new.cluster.ids[m1.group+1] = 'M1 macro'
new.cluster.ids[m2a.group+1] = 'M2a macro'
new.cluster.ids[m2c.group+1] = 'M2c macro'
new.cluster.ids[m0.group+1] = 'M0 macro'
# new.cluster.ids[fabp4.group+1] = 'Fabp4+ mono'
# new.cluster.ids[s100a8.group+1] = 'S100a8+ mono'
# new.cluster.ids[cd209a.group+1] = 'Cd209a+ mono'
# new.cluster.ids[xcr1.group+1] = 'Xcr1+ mono'
# new.cluster.ids[cd24a.group+1] = 'Cd24a+ mono'
# new.cluster.ids[slc8a1.group+1] = 'Slc8a1+ mono'
# new.cluster.ids[ccl5.group+1] = 'Ccl5+ mono'
# new.cluster.ids[ccl22.group+1] = 'Ccl22+ mono'
names(new.cluster.ids) <- levels(Idents(dt))

new.cluster.ids

dt <- RenameIdents(dt, new.cluster.ids)

options(repr.plot.height=6, repr.plot.width=6)
DimPlot(dt, reduction = 'umap', label = T, repel = T)

save(dt, file='./mono1229.Rdata')





m1.macro = subset(dt, idents = 'M1 macro')
m2.macro = subset(dt, idents = c('M2a macro', 'M2c macro'))
m0.macro = subset(dt, idents = 'M0 macro')
# all.macro = subset(dt, idents = c('M1 macro', 'M2 macro'))

dt1 = m2.macro@assays$RNA@data['Fabp4', ]
dt2 = m2.macro@assays$RNA@data['Slc27a1', ]
cor.test(dt1, dt2)
plot(dt1, dt2)



save(m1.macro, file='./sc_m1_from_mono.Rdata')
save(m2.macro, file='./sc_m2_from_mono.Rdata')
save(m0.macro, file='./sc_m0_from_mono.Rdata')

p1 = VlnPlot(m1.macro, features = 'M1.score1', pt.size = -1, sort = T, ncol=1, group.by = 'group') + NoLegend()
p2 = VlnPlot(m2.macro, features = 'M2.score1', pt.size = -1, sort = T, ncol=1, group.by = 'group') + NoLegend()
options(repr.plot.height=12, repr.plot.width=6)
plot_grid(plotlist = list(p1, p2), ncol = 1)

gg = p1$data
ggplot(gg, aes(x=ident, y=M1.score1)) + geom_boxplot()

options(repr.plot.height=8, repr.plot.width=12)
p1 = DotPlot(m1.macro, features = 'M1.score1', group.by = 'sample')
p1

gg = p1$data
ggg = gg[gg$features.plot == 'M1.score1', ]
ggg$group = "C"
ggg[grep('D', ggg$id), ]$group = 'D'
ggg[grep('GM', ggg$id), ]$group = 'GM'
ggg[grep('BBP', ggg$id), ]$group = 'BBP'
ggplot(ggg, aes(x=group, y=avg.exp.scaled)) + geom_boxplot()





genes = c('Slc44a1', 'Slc27a1', 'Il1b', 'Tnf', 'Fabp4', 'Apoa5', 'Pparg')
genes %in% rownames(data)

dt1 = m1.macro@assays$RNA@data['Slc44a1', m1.macro$group == 'GQD']
dt2 = m1.macro@assays$RNA@data['Slc44a1', m1.macro$group == 'D']

t.test(dt1, dt2)



options(repr.plot.height=12, repr.plot.width=12)
VlnPlot(m1.macro, features = genes, pt.size = 1, sort = T, group.by = 'group', ncol = 2)

options(repr.plot.height=12, repr.plot.width=12)
VlnPlot(m2.macro, features = genes, pt.size = 1, sort = T, group.by = 'group', ncol = 2)



options(repr.plot.height=8, repr.plot.width=12)
p1 = DotPlot(m1.macro, features = genes, group.by = 'sample')
p1

gg = p1$data
ggg = gg[gg$features.plot == 'Fabp4', ]
ggg$group = "C"
ggg[grep('D', ggg$id), ]$group = 'D'
ggg[grep('GM', ggg$id), ]$group = 'GM'
ggg[grep('BBP', ggg$id), ]$group = 'BBP'
ggplot(ggg, aes(x=group, y=avg.exp.scaled)) + geom_boxplot()



macro = subset(dt, idents=c('M1 macro', 'M0 macro', 'M2a macro', 'M2c macro'))

macro = FindVariableFeatures(macro, selection.method = "vst", nfeatures = 2000)
# data <- ScaleData(data, features = all.genes)
macro <- ScaleData(macro)

macro <- RunPCA(macro, verbose = FALSE, features = VariableFeatures(object = macro), npcs = 50)
macro <- FindNeighbors(macro, dims = 1:20, verbose = FALSE)
macro <- FindClusters(macro, verbose = FALSE, resolution = 1)
macro = RunUMAP(macro, dims = 1:20)

options(repr.plot.height=6, repr.plot.width=6)
DimPlot(macro, reduction = 'umap', label = T)

options(repr.plot.height=12, repr.plot.width=12)
FeaturePlot(macro, features = c('M1.score1', 'M2.score1', 'Mrc1', 'M0', 'M1', 'M2', 'nFeature_RNA', 'Tnf', 'Arg1'), 
            min.cutoff = 0, label = T, ncol = 3)

M1.like = 12
M1.group = 4
M2.group = c(3,8,9,11)

dtt = macro

new.cluster.ids = rep('M0 macro', length(levels(dtt)))
new.cluster.ids[M1.group+1] = 'M1 Tnf high'
new.cluster.ids[M2.group+1] = 'M2 macro'
new.cluster.ids[M1.like+1] = 'M1 Tnf low'
names(new.cluster.ids) <- levels(Idents(dtt))

new.cluster.ids

dtt <- RenameIdents(dtt, new.cluster.ids)

options(repr.plot.height=6, repr.plot.width=6)
DimPlot(dtt, reduction = 'umap', label = T, repel = T)

# options(repr.plot.height=10, repr.plot.width=8)
# DotPlot(dtt, features = c('M1.score1', 'M2.score1', genes), 
#         cols = colorRampPalette(c("blue", "red"))(4), dot.scale = 8, split.by = "group") + RotatedAxis()

options(repr.plot.height=10, repr.plot.width=8)
DotPlot(dtt, features = genes, 
        cols = colorRampPalette(c("blue", "red"))(4), dot.scale = 8, split.by = "group") + RotatedAxis()

table(dtt$group, Idents(dtt))

table(dtt$group)

(208+40)/3240
(329+89)/5013
(170+46)/4323
(111+39)/3304

(208+40)/2160 
(329+89)/2595
(170+46)/2490
(111+39)/1983

dtt$macro_type = Idents(dtt)

dtt@meta.data$macro_type2 = dtt@meta.data$macro_type
levels(dtt@meta.data$macro_type2) = c(levels(dtt@meta.data$macro_type2), 'M1 macro')
dtt@meta.data[dtt@meta.data$macro_type %in% c('M1 Tnf high', 'M1 Tnf low'), ]$macro_type2 = 'M1 macro'

dtt@meta.data$macro_group = paste0(dtt@meta.data$macro_type, '_', dtt@meta.data$group)
dtt@meta.data$macro_group2 = paste0(dtt@meta.data$macro_type2, '_', dtt@meta.data$group)

unique(dtt@meta.data$macro_group)
unique(dtt@meta.data$macro_group2)

save(dtt, file='./macro1229.Rdata')



dtt = macro

new.cluster.ids = rep('M0 macro', length(levels(dtt)))
new.cluster.ids[M1.group+1] = 'M1 Tnf high'
new.cluster.ids[M2.group+1] = 'M2 macro'
new.cluster.ids[M1.like+1] = 'M1 Tnf low'
names(new.cluster.ids) <- levels(Idents(dtt))

new.cluster.ids

dtt <- RenameIdents(dtt, new.cluster.ids)

m1.from_macro = subset(dtt, idents = c('M1 Tnf high', 'M1 Tnf low'))
m2.from_macro = subset(dtt, idents = 'M2 macro')
m0.from_macro = subset(dtt, idents = 'M0 macro')

save(m1.from_macro, file='./sc_m1_from_macro.Rdata')
save(m2.from_macro, file='./sc_m2_from_macro.Rdata')
save(m0.from_macro, file='./sc_m0_from_macro.Rdata')







comparison1 = data.frame(front=c('M1 Tnf high_D', 'M1 Tnf high_GQD', 'M1 Tnf high_BBP', 'M1 Tnf high_GQD'), 
                         back=c('M1 Tnf high_C', 'M1 Tnf high_D', 'M1 Tnf high_D', 'M1 Tnf high_BBP'))

comparison2 = data.frame(front=c('M1 Tnf low_D', 'M1 Tnf low_GQD', 'M1 Tnf low_BBP', 'M1 Tnf low_GQD'), 
                         back=c('M1 Tnf low_C', 'M1 Tnf low_D', 'M1 Tnf low_D', 'M1 Tnf low_BBP'))

comparison3 = data.frame(front=c('M1 Tnf high_C', 'M1 Tnf high_D', 'M1 Tnf high_BBP', 'M1 Tnf high_GQD'), 
                         back=c('M1 Tnf low_C', 'M1 Tnf low_D', 'M1 Tnf low_BBP', 'M1 Tnf low_GQD'))

Idents(dtt) = dtt$macro_group

diff.list1 = list()
for(i in 1:dim(comparison1)[1]){
    diff = FindMarkers(dtt, `ident.1` = comparison1[i,1], `ident.2` = comparison1[i,2])
    diff.list1[[i]] = diff
    names(diff.list1)[i] = paste0(comparison1[i,1], '_', comparison1[i,2])
}

diff.list2 = list()
for(i in 1:dim(comparison2)[1]){
    diff = FindMarkers(dtt, `ident.1` = comparison2[i,1], `ident.2` = comparison2[i,2])
    diff.list2[[i]] = diff
    names(diff.list2)[i] = paste0(comparison2[i,1], '_', comparison2[i,2])
}

diff.list3 = list()
for(i in 1:dim(comparison3)[1]){
    diff = FindMarkers(dtt, `ident.1` = comparison3[i,1], `ident.2` = comparison3[i,2])
    diff.list3[[i]] = diff
    names(diff.list3)[i] = paste0(comparison3[i,1], '_', comparison3[i,2])
}

diff.list = list(high=diff.list1, low=diff.list2, h_vs_l=diff.list3)

save(diff.list, file='./test_m1_h_l.Rdata')





comparison1 = data.frame(front=c('M1 macro_D', 'M1 macro_GQD', 'M1 macro_BBP', 'M1 macro_GQD'), 
                         back=c('M1 macro_C', 'M1 macro_D', 'M1 macro_D', 'M1 macro_BBP'))

Idents(dtt) = dtt$macro_group2

diff.list1 = list()
for(i in 1:dim(comparison1)[1]){
    diff = FindMarkers(dtt, `ident.1` = comparison1[i,1], `ident.2` = comparison1[i,2])
    diff.list1[[i]] = diff
    names(diff.list1)[i] = paste0(comparison1[i,1], '_', comparison1[i,2])
}

diff.list = diff.list1

save(diff.list, file='./test_m1.Rdata')



options(repr.plot.height=20, repr.plot.width=16)a
VlnPlot(dtt, features = genes, pt.size = 1, sort = T, ncol = 1, split.by = 'group')





markers = FindAllMarkers(dt, only.pos = T)

markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

DoHeatmap(dt, features = top10$gene) + NoLegend()





options(repr.plot.height=10, repr.plot.width=12)
VlnPlot(m0.macro, features = genes, pt.size = 1, sort = T, ncol = 2, group.by = 'group')

t1 = m0.macro@assays$RNA@data['Slc44a1', m0.macro$group == 'BBP']
t2 = m0.macro@assays$RNA@data['Slc44a1', m0.macro$group == 'D']
t3 = m0.macro@assays$RNA@data['Slc44a1', m0.macro$group == 'GQD']
t.test(t1, t2)
t.test(t3, t2)
length(t1[(t1 > 0)])/length(t1)
length(t2[(t2 > 0)])/length(t2)
length(t3[(t3 > 0)])/length(t3)





dt1 = m1.macro@assays$RNA@data['Tnf', ]
dt2 = m1.macro@assays$RNA@data['Lpl', ]
# dt1 = m1.macro@assays$RNA@data['Tnf', ]
# dt2 = m1.macro@assays$RNA@data['Lpl', ]
# cor.test(dt1, dt2)
dt2.high = names(dt2[dt2 > 0])
dt1.high = names(dt1[dt1 > 0])
both.high = intersect(dt1.high, dt2.high)
cor.test(dt1[both.high], dt2[both.high])
plot(dt1[both.high], dt2[both.high])




