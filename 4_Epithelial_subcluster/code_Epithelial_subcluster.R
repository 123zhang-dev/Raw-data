library(Seurat)
library(tidyverse)
library(sctransform)
library(harmony)
library(glmGamPoi)
#########################################################
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
Epithelial <- CellCycleScoring(Epithelial, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
Epithelial = SCTransform(Epithelial, method="glmGamPoi", vars.to.regress=c("nCount_RNA", "percent.mito", "S.Score", "G2M.Score"), verbose=FALSE)
Epithelial = RunPCA(Epithelial, verbose=FALSE)
Epithelial = RunHarmony(Epithelial, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(Epithelial)
Epithelial <- FindNeighbors(Epithelial, dims=1:30, reduction="harmony")
Epithelial <- RunUMAP(Epithelial, dims=1:30, reduction="harmony")
p = DimPlot(Epithelial, reduction="umap", group.by="Sample", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("UMAP_Sample.pdf", p, height=6, width=6)
######################################################### 划分所有细胞的亚群时，resolution=0.1
mydata <- FindClusters(Epithelial, resolution=0.1)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors, label.size=5)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "./All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
############## 用小提琴图来检测marker基因在细胞亚群之间的分布
VlnPlot(mydata, features=c("S100A9"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())
# 0: JSRP1, HPN, CDH16, NKAIN4, WNT7A, RASIP1, ADGRE2, NEURL3, PIK3AP1, CBX6, APOBEC3B
# 1: BEX3, FOLR1, LTF, EYA2, SLC39A8, CHGA, LDHB, BPIFB1, STC1, PKHD1L1
# 2: BTNL3, PPP1R1B, PRAC1, SULT1B1, MOGAT3, SULT1E1, ALDOB, CDX1, BTNL8, CHP2, HMGCS2, MEP1A
colors = c("#FB8072", "#80B1D3", "#FDB462")
cell_label = c("Epithelial cells 1", "Epithelial cells 2", "Epithelial cells 3")
#################################################################################
## 给细胞的标签命名
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
saveRDS(mydata, "./mydata_cluster.rds")
UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)+theme(text=element_text(family="Times"))
genes = c("JSRP1", "HPN", "CDH16", "NKAIN4", "WNT7A", "RASIP1", "ADGRE2", "NEURL3", "BEX3", "FOLR1", "LTF", "EYA2", "SLC39A8", "CHGA", "LDHB", "BPIFB1", "BTNL3", "PPP1R1B", "PRAC1", "SULT1B1", "MOGAT3", "SULT1E1", "ALDOB", "CDX1")
p = DotPlot(mydata, features=genes)+coord_flip()+theme_minimal()+scale_color_gradientn(colors=c("dodgerblue", "white", "hotpink"))+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=13, face="bold"), axis.text.y=element_text(face="bold", size=14), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("marker_dotplot.pdf", p, width=7, height=6)

# marker基因的UMAP映射图
set = c("FOXP3", "GNLY", "TFF3", "CEACAM6", "SPP1", "SFRP4", "IGHG4", "BANK1")
FeaturePlot(mydata, features=set, cols=c("snow", "hotpink"), ncol=4)

set = c("IL7R", "GNLY", "MZB1", "LYZ", "MS4A1", "COL1A1", "VWF", "TPSAB1")


#####################################################################  配对柱状图
Type_label = c("Control", "24h Treated", "72h Treated")
bar$Type = factor(bar$Type, levels=Type_label)
bar = bar %>% group_by(Type) %>% mutate(percent=100*n/sum(n))

ggplot(data=bar2, aes(x=cell_type, y=percent, fill=Type))+
geom_bar(stat="identity", position=position_dodge())+
scale_fill_manual(values=c("burlywood1", "violet"))+theme_minimal()+geom_text(aes(label=percent), vjust=-0.2)+
theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=13, face="bold"), axis.text.y=element_text(face="bold", size=14), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
########################## 画配对秩和检验的箱线图
library(ggplot2)
library(ggpubr)

df = read.table("./celltype_number_percent.txt", header=T, sep="\t")
ggplot(df, aes(x=reorder(cell_type, -percent, sum), y=percent, fill=Type))+
scale_fill_manual(values=c("orange1", "greenyellow"))+
geom_boxplot(outlier.size=0.1, width=0.3)+
theme_bw()+
stat_compare_means(aes(group=Type), label="p.signif", method="t.test")+
theme(axis.text.x=element_text(angle=15, hjust=1, face="bold", size=10), axis.text.y=element_text(face="bold", size=10), axis.title.x=element_blank())


	



