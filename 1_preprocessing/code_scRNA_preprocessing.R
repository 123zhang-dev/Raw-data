library(Seurat)
library(tidyverse)
library(sctransform)
library(harmony)
library(glmGamPoi)
###############################################################################
assays <- dir("./0_rawdata/")
dir <- paste0("./0_rawdata/", assays)
# 
samples_name = assays
# 
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-") 
  }
}
### 
names(scRNAlist) <- samples_name
# 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<7000&percent.mito<10)  #
VlnPlot(scRNA, features=c("nFeature_RNA"), pt.size=0, cols=colors)
#########################################################################################################
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
scRNA <- CellCycleScoring(scRNA, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
scRNA = SCTransform(scRNA, method="glmGamPoi", vars.to.regress=c("nCount_RNA", "percent.mito", "S.Score", "G2M.Score"), verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
scRNA = RunHarmony(scRNA, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims=1:30, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:30, reduction="harmony")
p = DimPlot(scRNA, reduction="umap", group.by="Sample", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("UMAP_Sample.pdf", p, height=6, width=6)
#########################################################
mydata <- FindClusters(scRNA, resolution=0.1)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors, label.size=5)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "./All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
############## 
VlnPlot(mydata, features=c("S100A9"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())
# 0:   Regulatory T cells: FOXP3, CTLA4, ICOS, CCR7, SELL
# 1:   Cytotoxic NK/T cells: GNLY, NKG7, GZMA, CD8A, CD8B
# 2:   Epithelial cells: TFF3, KRT8, AGR2, GDF15, CLDN4, KRT17, CEACAM6
# 3:   Myeloid cells: SPP1, APOC1, FCER1G, S100A8, G0S2
# 4:   Myofibroblasts: SFRP4, COL1A1, LUM, TAGLN, ACTA2
# 5,6: Plasma B cells: IGHG4, IGLL5, PNOC
# 7:   B cells: MS4A1, BANK1, CD19, FCRLA
cell_label = c("Regulatory T cells", "Cytotoxic NK/T cells", "Epithelial cells", "Myeloid cells",
"Myofibroblasts", "Plasma B cells", "Plasma B cells", "B cells")
#################################################################################
## 
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
colors = c4a("tableau.10", 9)
saveRDS(mydata, "./mydata_cluster.rds")
colors = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977", "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A", "#39737C", "#86B4A9", "#82853B", "#CCC94D")
UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)+theme(text=element_text(family="Times"))
genes = c("FOXP3", "CTLA4", "ICOS", "GNLY", "NKG7", "CD8A", "KRT8", "AGR2", "GDF15", "CLDN4", "CEACAM6", "SPP1", "APOC1", "S100A8", "G0S2", "SFRP4", "COL1A1", "TAGLN", "IGHG4", "IGLL5", "MS4A1", "BANK1", "CD19")
p = DotPlot(mydata, features=genes)+coord_flip()+theme_minimal()+scale_color_gradientn(colors=c("dodgerblue", "white", "hotpink"))+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=13, face="bold"), axis.text.y=element_text(face="bold", size=14), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("marker_dotplot.pdf", p, width=7, height=6)

# marker
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
########################## 
library(ggplot2)
library(ggpubr)

df = read.table("./celltype_number_percent.txt", header=T, sep="\t")
ggplot(df, aes(x=reorder(cell_type, -percent, sum), y=percent, fill=Type))+
scale_fill_manual(values=c("orange1", "greenyellow"))+
geom_boxplot(outlier.size=0.1, width=0.3)+
theme_bw()+
stat_compare_means(aes(group=Type), label="p.signif", method="t.test")+
theme(axis.text.x=element_text(angle=15, hjust=1, face="bold", size=10), axis.text.y=element_text(face="bold", size=10), axis.title.x=element_blank())


	



