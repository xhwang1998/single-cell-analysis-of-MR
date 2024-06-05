library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(stringr)
library(data.table)
library(RColorBrewer)
library(readr)
library(Cairo)
library(clustree)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)
library(ggpubr)

setwd("./Analysis/week20220620/data/Mye")
Mye.data <- subset(useful.cells,Annotation=="Myeloid")
Mye.data <- fuwai.subcluster.integration(data=Mye.data,prefix="Mye")#nPCs=12
fuwai.subcluster.plot(data=Mye.data,prefix="Mye")
DimPlot(object = Mye.data,
        reduction = "umap",
        pt.size = 0.1,
        label = TRUE,
        #label = FALSE,
        raster=FALSE,
        label.size = 6,
        group.by= "integrated_snn_res.0.2"
) + ggtitle(label = 'Umap_Annotation') + coord_fixed(1:1)

#############################################################
#find.markers
DefaultAssay(Mye.data) <- "RNA"
top5.Mye <- list()
cellcount.list <- list()
for (res in c(seq(0.1,0.5,0.05))){
  top5.res <- find.markers(data=Mye.data,res=paste0("integrated_snn_res.",res),cell.name=paste0("Mye_res.",res),n=5)
  top5.Mye[[paste0("Res.",res)]] <- top5.res
  cellcount.list[[paste0("Res.",res)]] <- table(Mye.data$Sample,Mye.data@meta.data[[paste0("integrated_snn_res.",res)]])
}
write.xlsx(cellcount.list,file="Mye_ClusterCount.xlsx",rowNames =TRUE)
##########################################################
#Annotation
#Mye老师确定分辨率选择0.15
Idents(Mye.data) <- "integrated_snn_res.0.15"
Mye.data <- RenameIdents(Mye.data,
                         '0'='Mye0',
                         '1'='Mye1',
                         '2'='Mye2',
                         '3'='Mye3',
                         '4'='Mye4',
                         '5'='Mye5')
Mye.data$sub.anno <- Mye.data@active.ident
saveRDS(Mye.data,"Mye_Integration_UMAP_annotation.rds")

#####################################################################
#GO-BP and KEGG analysis
setwd('./Analysis/week20220704/data/Mye')
mye.deg.raw <- read.csv("Mye_res.0.15_Find_AllMarkers.tab",sep="\t")
head(mye.deg.raw,3)

mye.deg.raw %>% filter(avg_logFC >0.5) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> Mye.top50
Mye.data <- ScaleData(Mye.data)
pdf('Mye_subcluster_heatmapDEGs.pdf')
DoHeatmap(Mye.data, features = Mye.top50$gene,group.by = "sub.anno",label=FALSE)+theme(axis.text.y = element_text(size = 2))
dev.off()

##########################################
#enrichment analysis
mye.deg.raw %>% filter(pct.1 >0.25) %>% filter(p_val < 0.05) -> mye.deg.filtered
anno.mye <- c("Mye0","Mye1","Mye2","Mye3","Mye4","Mye5")
names(anno.mye) <- c("0","1","2","3","4","5")
mye.deg.filtered$anno <- anno.mye[as.character(mye.deg.filtered$cluster)]
Mye.enrich <- enrich_fuwai_new(degs=mye.deg.filtered,OrgDb="org.Hs.eg.db",organism="hsa",prefix="Mye")
######################################################################
#Group compare
library(ggpubr)
Mye.group <- as.data.frame(prop.table(table(Mye.data$Sample,Mye.data$sub.anno),margin = 1))
Mye.group$Group <- ifelse(Mye.group$Var1 %in% c("MR1","MR2","MR3"),"MR","NC")
Mye.group$ratio <- round(Mye.group$Freq*100,1)
ggpubr::compare_means(ratio~Group,Mye.group,
                      method = "wilcox.test",paired = FALSE,
                      group.by = "Var2")
######################################################################
#SCENIC
Mye.exprMat <- t(as.matrix(Mye.data@assays$RNA@counts))
write.table(Mye.exprMat,"Mye_RawCounts.tsv",sep="\t",quote=F)

#pySCENIC运行完后绘图
setwd("./Analysis/week20220704/data/Mye/scenic")
write.table(data.frame(Mye.data$sub.anno),"Mye.celltype",sep="\t",quote=F)
scenic.plot(adata=Mye.data,prefix="Mye")

#############################################################################
#cell cycle score
#Seurat CellCycleScoring
setwd("./Analysis/week20220704/data/Mye/cecllcycle")
s.genes <- CaseMatch(search=cc.genes$s.genes ,match=rownames(Mye.data))
g2m.genes <- CaseMatch(search=cc.genes$g2m.genes,match=rownames(Mye.data))
set.seed(5)
Mye.cyclescore <- CellCycleScoring(Mye.data, s.features = s.genes, g2m.features = g2m.genes)
saveRDS(Mye.cyclescore,'Mye_Umap_anno_CellCycleScore.rds')
score.data <- data.frame(Mye.cyclescore$sub.anno,Mye.cyclescore$S.Score,Mye.cyclescore$G2M.Score,Mye.cyclescore$Phase)
colnames(score.data) <- c("Cluster","S.Score","G2M.Score","Phase")
write.xlsx(score.data,'Mye_CellCycleScores.xlsx',rowNames =TRUE)

#scran包的cyclone函数计算
assignments <- readRDS('./Mye_cyclone_CellCycle.rds')
anno.df <- Mye.data@meta.data[,c('Sample','Group','sub.anno')]
anno.df$Phases <- assignments$phases
mye.scores <- data.frame(t(assignments$scores))
colnames(mye.scores) <- rownames(anno.df)

#https://sunpma.com/other/rgb/
#p1 <- DimPlot(VEC.data)
#unique(factor(ggplot_build(p1)$data[[1]]$colour))
ann_colors = list(
  Phases = c(G1="#0CDA5F",S="#8BB4ED",G2M="#FC9489"),
  Group = c('MR'="#7FC97F",'NC'="#BEAED4"),
  Sample = c('MR1'="#8DD3C7",
             'MR2'="#FFFFB3",
             'MR3'="#BEBADA",
             'NC1'="#FB8072",
             'NC2'="#80B1D3",
             'NC3'="#FDB462"),
  sub.anno = c('Mye0'='#F8766D',
               'Mye1'='#B79F00',
               'Mye2'='#00BA38',
               'Mye3'='#00BFC4',
               'Mye4'='#619CFF',
               'Mye5'='#F564E3'
               )
)
pheatmap(mye.scores,
         show_colnames = F,
         #cluster_cols = FALSE, 
         cluster_rows = FALSE,
         annotation_col= anno.df,
         annotation_colors = ann_colors,
         filename="Mye_CellScore_pheatmap.pdf",
         width=15,height = 8,
         cellwidth = 0.2, cellheight = 130)

############################################################################################################
#Group 差异分析
setwd('./Analysis/week20220718/data/Mye')
Idents(seurat.data) <- 'Group'
markers.raw <- FindAllMarkers(seurat.data,test.use = "t", min.cells.group = 0, min.pct=0, logfc.threshold = -Inf)
mye.markers.raw <- readRDS('Mye.Group_Allmarkers.rds')
fuwai.deg.enrich(markers.raw = mye.markers.raw, prefix = 'Mye')
valcano.plot(markers=mye.markers.raw,prefix='Mye')
