library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(stringr)
library(data.table)
library(RColorBrewer)
library(readr)
library(angrycell,lib.loc='/BGFS1/projectdata/jasper/database/Rlib/R_4.1.2')
library(Cairo)
library(clustree)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)Lym
library(ggpubr)

setwd("./Analysis/week20220620/data/Lym")
Lym.data <- subset(useful.cells,Annotation=="Lymphocyte")
Lym.data <- fuwai.subcluster.integration(data=Lym.data,prefix="Lym")#nPCs=9
Lym.data <- FindClusters(object = Lym.data,
                         verbose = F,
                         resolution = c(seq(0.6,0.9,0.1)))
clustree(Lym.data@meta.data, prefix = "integrated_snn_res.")
fuwai.subcluster.plot(data=Lym.data,prefix="Lym")

DimPlot(object = Lym.data,
        reduction = "umap",
        pt.size = 0.1,
        label = TRUE,
        #label = FALSE,
        raster=FALSE,
        label.size = 6,
        group.by= "integrated_snn_res.0.9"
) + ggtitle(label = 'Umap_Annotation') + coord_fixed(1:1)

#############################################################
#find.markers
DefaultAssay(Lym.data) <- "RNA"
top5.Lym <- list()
cellcount.list <- list()
for (res in c(seq(0.1,0.5,0.05))){
  top5.res <- find.markers(data=Lym.data,res=paste0("integrated_snn_res.",res),cell.name=paste0("Lym_res.",res),n=5)
  top5.Lym[[paste0("Res.",res)]] <- top5.res
  cellcount.list[[paste0("Res.",res)]] <- table(Lym.data$Sample,Lym.data@meta.data[[paste0("integrated_snn_res.",res)]])
}
write.xlsx(cellcount.list,file="Lym_ClusterCount.xlsx",rowNames =TRUE)

#Classical markers
markers.list <- c("MS4A1","CD79A","CD79B","CD3D","CD3E","CD3G","CD4","CD8A","CD8B","NCR1","NCAM1", "NKG7", "KLRD1","KLRF1","LUM","COL1A1","COL3A1","SOX9","ACTA2","MYH11")
pdf("Lym_classical_markers.pdf",h=9,w=15)
FeaturePlot(Lym.data,features = markers.list , order=TRUE,ncol=5,pt.size=0.1)
dev.off()
DotPlot(Lym.data,features = c("MS4A1","CD79A","CD79B","CD3D","CD3E","CD3G","CD4","IL7R","CD8A","CD8B","NCR1","NCAM1", "NKG7", "KLRD1","KLRF1"))

Lym.data <- readRDS('Lym_Data-Integration_UMAP.rds')
pdf('Lym_VlnPlot_classical_markers.pdf',h=2.5,w=9.5)
for (res in c(seq(0.1,0.5,0.05))){
  cluster.n <- length(levels(Lym.data@meta.data[[paste0("integrated_snn_res.",res)]]))
  Idents(Lym.data) <- paste0("integrated_snn_res.",res)
  pi <- vioplot2(Lym.data,paths=markers.list,facet.toward = "col",remove.axis.x = T,panel.spacing = 0)+
    theme(plot.margin = margin(t = 40-cluster.n*3, b = 40-cluster.n*3, r=20, l=20),
          plot.title = element_text(hjust = 0.5))+
    labs(x="",y="",title=paste0("Resolution: ",res))
  print(pi)
}
dev.off()

##########################################################
#Annotation
#Lym老师确定分辨率选择0.3
Idents(Lym.data) <- "integrated_snn_res.0.3"
Lym.data <- RenameIdents(Lym.data,
                         '0'='Lym0',
                         '1'='Lym1',
                         '2'='Lym2',
                         '3'='Lym3',
                         '4'='Lym4',
                         '5'='Lym5')
Lym.data$sub.anno <- Lym.data@active.ident
saveRDS(Lym.data,"Lym_Integration_UMAP_annotation.rds")

pdf("Lym_res.0.3_anno.cluster.pdf",w=10,h=6)
DimPlot(object = Lym.data,
        reduction = "umap",
        pt.size = 0.5,
        label = TRUE,
        #label = FALSE,
        raster=FALSE,
        label.size = 5,
        #group.by= "Group"
) + ggtitle(label = 'Umap_Annotation') + coord_fixed(1:1)
dev.off()

#####################################################################
#GO-BP and KEGG analysis
setwd('./Analysis/week20220704/data/Lym')
lym.deg.raw <- read.csv("Lym_res.0.3_Find_AllMarkers.tab",sep="\t")
head(lym.deg.raw,3)

lym.deg.raw %>% filter(avg_logFC >0.5) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> Lym.top50
Lym.data <- ScaleData(Lym.data)
pdf('Lym_subcluster_heatmapDEGs.pdf')
DoHeatmap(Lym.data, features = Lym.top50$gene,group.by = "sub.anno",label=FALSE)+theme(axis.text.y = element_text(size = 2))
dev.off()

##########################################
#enrichment analysis
lym.deg.raw %>% filter(pct.1 >0.25) %>% filter(p_val < 0.05) -> lym.deg.filtered
anno.lym <- c("Lym0","Lym1","Lym2","Lym3","Lym4","Lym5")
names(anno.lym) <- c("0","1","2","3","4","5")
lym.deg.filtered$anno <- anno.lym[as.character(lym.deg.filtered$cluster)]
Lym.enrich <- enrich_fuwai_new(degs=lym.deg.filtered,OrgDb="org.Hs.eg.db",organism="hsa",prefix="Lym")

######################################################################
#Group compare
library(ggpubr)
Lym.group <- as.data.frame(prop.table(table(Lym.data$Sample,Lym.data$sub.anno),margin = 1))
Lym.group$Group <- ifelse(Lym.group$Var1 %in% c("MR1","MR2","MR3"),"MR","NC")
Lym.group$ratio <- round(Lym.group$Freq*100,1)
ggpubr::compare_means(ratio~Group,Lym.group,
                      method = "wilcox.test",paired = FALSE,
                      group.by = "Var2")
######################################################################
#SCENIC
Lym.exprMat <- t(as.matrix(Lym.data@assays$RNA@counts))
write.table(Lym.exprMat,"Lym_RawCounts.tsv",sep="\t",quote=F)

#pySCENIC运行完后绘图
setwd("./Analysis/week20220704/data/Lym/scenic")
write.table(data.frame(Lym.data$sub.anno),"Lym.celltype",sep="\t",quote=F)
scenic.plot(adata=Lym.data,prefix="Lym")

############################################################################################################
#Group 差异分析
Idents(seurat.data) <- 'Group'
markers.raw <- FindAllMarkers(seurat.data,test.use = "t", min.cells.group = 0, min.pct=0, logfc.threshold = -Inf)
setwd('./Analysis/week20220718/data/Lym')
lym.markers.raw <- readRDS('Lym.Group_Allmarkers.rds')
fuwai.deg.enrich(markers.raw = lym.markers.raw, prefix = 'Lym')
valcano.plot(markers=mye.markers.raw,prefix='Mye')
