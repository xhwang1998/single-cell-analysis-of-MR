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

setwd("./Analysis/week20220620/data/VEC")
VEC.data <- subset(useful.cells,Annotation=="VEC")
sub.data <- VEC.data
prefix <- "VEC"
DefaultAssay(sub.data) <- "RNA"
data.list <- SplitObject(sub.data,split.by = "Sample")
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = data.list)
#细胞数太少时FindIntegrationAnchors和IntegrateData会报错
group.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features,k.filter = NA)
sub.data <- IntegrateData(anchorset = group.anchors, k.weight = 30)
DefaultAssay(sub.data) <- "integrated"
sub.data <- ScaleData(object = sub.data,vars.to.regress = "percent.mt",verbose = FALSE)
sub.data <- RunPCA(object = sub.data, verbose = F)
P_PCA <- DimPlot(object = sub.data, reduction = "pca",group.by = "Sample") + coord_fixed(1:1)
P_Elbow <- ElbowPlot(sub.data)
pdf(file=paste0(prefix,"_DimensionalReduction_PCA.pdf"),width = 16,height = 9)
print(P_PCA)
print(P_Elbow)
dev.off()

x <- sub.data
pct<-x[["pca"]]@stdev/sum(x[["pca"]]@stdev)*100
cumu<-cumsum(pct)
co1<-which(cumu >80 & pct<5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
nPC <- min(co1,co2)
print(paste0("PCA & nPCs Selected Finished. ",nPC," PCs Selected"))
rm(x,co1,co2,cumu,pct)
print(nPC)
sub.data <- FindNeighbors(object = sub.data,
                          reduction = "pca", 
                          dims = 1:nPC,
                          verbose = F)
sub.data <- FindClusters(object = sub.data,
                         verbose = F,
                         resolution = c(seq(0.1,0.5,0.05)))
p_clust <- clustree(sub.data@meta.data, prefix = "integrated_snn_res.")
pdf(file=paste0(prefix,'_clustree_plot.pdf'),width = 16,height = 9)
print(p_clust)
dev.off()
sub.data <- RunUMAP(object = sub.data, 
                    dims = 1:nPC)
saveRDS(sub.data,file = paste0(prefix,"_Data-Integration_UMAP.rds"))
###########################################################
nPC=9
VEC.data <- sub.data
DimPlot(object = VEC.data,
        reduction = "umap",
        pt.size = 0.00001,
        #label = TRUE,
        label = FALSE,
        raster=FALSE,
        group.by= "Group"
) + ggtitle(label = 'Umap_Group') + coord_fixed(1:1)
fuwai.subcluster.plot(data=VEC.data,prefix="VEC.cells")
#############################################################
#find.markers
DefaultAssay(VEC.data) <- "RNA"
top5.VEC <- list()
cellcount.list <- list()
for (res in c(seq(0.1,0.5,0.05))){
  top5.res <- find.markers(data=VEC.data,res=paste0("integrated_snn_res.",res),cell.name=paste0("VEC_res.",res),n=5)
  top5.VEC[[paste0("Res.",res)]] <- top5.res
  cellcount.list[[paste0("Res.",res)]] <- table(VEC.data$Sample,VEC.data@meta.data[[paste0("integrated_snn_res.",res)]])
}
write.xlsx(cellcount.list,file="VEC_ClusterCount.xlsx",rowNames =TRUE)
##########################################################
#Annotation
#VEC老师确定分辨率选择0.2
Idents(VEC.data) <- "integrated_snn_res.0.2"
VEC.data <- RenameIdents(VEC.data,
                         '0'="VEC0",
                         '1'="VEC1",
                         '2'="VEC2")
VEC.data$sub.anno <- VEC.data@active.ident
saveRDS(VEC.data,"VEC_Integration_UMAP_annotation.rds")
pdf("VEC_res.0.2_anno.cluster.pdf",w=10,h=6)
DimPlot(object = VEC.data,
        reduction = "umap",
        pt.size = 1,
        label = TRUE,
        #label = FALSE,
        raster=FALSE,
        label.size = 6,
        #group.by= "Group"
) + ggtitle(label = 'Umap_Annotation') + coord_fixed(1:1)
dev.off()

#####################################################################
#GO-BP and KEGG analysis
setwd('.//Analysis/week20220704/data/VEC')
vec.deg.raw <- read.csv("VEC_res.0.2_Find_AllMarkers.tab",sep="\t")
head(vec.deg.raw,3)

vec.deg.raw %>% filter(avg_logFC >0.5) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> VEC.top50
VEC.data <- ScaleData(VEC.data)
pdf('VEC_subcluster_heatmapDEGs.pdf')
DoHeatmap(VEC.data, features = VEC.top50$gene,group.by = "sub.anno",label=FALSE)+theme(axis.text.y = element_text(size = 2))
dev.off()

##########################################
#enrichment analysis
vec.deg.raw %>% filter(pct.1 >0.25) %>% filter(p_val < 0.05) -> vec.deg.filtered
anno.vec <- c("VEC0","VEC1","VEC2")
names(anno.vec) <- c("0","1","2")
vec.deg.filtered$anno <- anno.vec[as.character(vec.deg.filtered$cluster)]
VEC.enrich <- enrich_fuwai_new(degs=vec.deg.filtered,OrgDb="org.Hs.eg.db",organism="hsa",prefix="VEC")
for (logfc in c(seq(0.25,1,0.25))){
  flag <- paste0("logfc>",logfc)
  flag <- "logfc>1"
  p1 <- clusterProfiler::dotplot(VEC.enrich$go[[flag]], showCategory=5, includeAll=FALSE,label_format=100)+
    labs(x='',y='',title="GO Enrichment Analysis")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",colour = "black",size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.margin = margin(t = 20, b = 20, r = 20,  l = 20))
  ggsave(paste0("VEC_GOBP_",flag,"_top5Dotplot.pdf"),p1,width = 7.3,height = 5.5)
  
  p2 <- clusterProfiler::dotplot(VEC.enrich$kegg[[flag]], showCategory=5, includeAll=FALSE,label_format=100)+
    labs(x='',y='',title="KEGG Enrichment Analysis")+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",colour = "black",size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.margin = margin(t = 10, b = 10, r = 80,  l = 80))
  ggsave(paste0("VEC_KEGG_",flag,"_top5Dotplot.pdf"),p2,width = 8.2,height = 4.5)
}

######################################################################
#Group compare
library(ggpubr)
VEC.group <- as.data.frame(prop.table(table(VEC.data$Sample,VEC.data$sub.anno),margin = 1))
VEC.group$Group <- ifelse(VEC.group$Var1 %in% c("MR1","MR2","MR3"),"MR","NC")
VEC.group$ratio <- round(VEC.group$Freq*100,1)
######################################################################
#SCENIC
VEC.exprMat <- t(as.matrix(VEC.data@assays$RNA@counts))
write.table(VEC.exprMat,"VEC_RawCounts.tsv",sep="\t",quote=F)

#pySCENIC运行完后绘图
setwd(".//Analysis/week20220704/data/VEC/scenic")
write.table(data.frame(VEC.data$sub.anno),"VEC.celltype",sep="\t",quote=F)
#/BGFS1/home/yanyl/miniconda3/envs/pyscenic/bin/python ../../../script/AUC_CellToCluster.py --in_dir . --out_dir . --prefix VEC
scenic.plot(adata=VEC.data,prefix="VEC")

#############################################################################
#cell cycle score
#Seurat CellCycleScoring
setwd(".//Analysis/week20220704/data/VEC/cecllcycle")
s.genes <- CaseMatch(search=cc.genes$s.genes ,match=rownames(VEC.data))
g2m.genes <- CaseMatch(search=cc.genes$g2m.genes,match=rownames(VEC.data))
set.seed(5)
VEC.cyclescore <- CellCycleScoring(VEC.data, s.features = s.genes, g2m.features = g2m.genes)
saveRDS(VEC.cyclescore,'VEC_Umap_anno_CellCycleScore.rds')
score.data <- data.frame(VEC.cyclescore$sub.anno,VEC.cyclescore$S.Score,VEC.cyclescore$G2M.Score,VEC.cyclescore$Phase)
colnames(score.data) <- c("Cluster","S.Score","G2M.Score","Phase")
write.xlsx(score.data,'VEC_CellCycleScores.xlsx',rowNames =TRUE)

#scran包的cyclone函数计算
assignments <- readRDS('./VEC_cyclone_CellCycle.rds')
anno.df <- VEC.data@meta.data[,c('Sample','Group','sub.anno')]
anno.df$Phases <- assignments$phases
vec.scores <- data.frame(t(assignments$scores))
colnames(vec.scores) <- rownames(anno.df)

#p1 <- DimPlot(V.data,group.by = 'Sample')
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
  sub.anno = c('VEC0'='#F8766D',
               'VEC1'='#00BA38',
               'VEC2'='#619CFF')
)
pheatmap(vec.scores,
         show_colnames = F,
         #cluster_cols = FALSE, 
         cluster_rows = FALSE,
         annotation_col= anno.df,
         annotation_colors = ann_colors,
         filename="VEC_CellScore_pheatmap.pdf",
         width=15,height = 8,
         cellwidth = 1, cellheight = 130)

################################################################################
#Group 差异分析
setwd('.//Analysis/week20220718/data/VEC')
Idents(seurat.data) <- 'Group'
markers.raw <- FindAllMarkers(seurat.data,test.use = "t", min.cells.group = 0, min.pct=0, logfc.threshold = -Inf)
vec.markers.raw <- readRDS('VEC.Group_Allmarkers.rds')
fuwai.deg.enrich(markers.raw = vec.markers.raw, prefix = 'VEC')
valcano.plot(markers=vec.markers.raw,prefix='VEC')

