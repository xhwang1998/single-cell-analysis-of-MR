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
library(limma)
devtools::install_github("aertslab/SCENIC") 
library("AUCell",lib.loc = "/BGFS1/home/yanyl/miniconda3/envs/pyscenic/lib/R/library")
#install.packages("devtools")
library(SCENIC)
library("ComplexHeatmap")
                         
setwd("./Analysis/week20220620/data/VIC")
#useful.cells <- readRDS('./Analysis/week20220613/data/Rdata/Integration.Umap.annotation_discardLowQuality.rds')
colors.sample <- brewer.pal(9,"Set3")[1:6]
colors.group <- brewer.pal(8,"Accent")[1:2]
VIC.data <- subset(useful.cells,Annotation=="VIC")
VIC.data <- fuwai.subcluster.integration(data=VIC.data,prefix="VIC.Cells")#nPCs=15

DefaultAssay(VIC.data) <- "RNA"
top5.list <- list()
cellcount.list <- list()
for (res in c(seq(0.1,0.5,0.05))){
  top5.res <- find.markers(data=VIC.data,res=paste0("integrated_snn_res.",res),cell.name=paste0("VIC_res.",res),n=5)
  top5.list[[paste0("Res.",res)]] <- top5.res
  dot.list[[paste0("Res.",res)]] <- dot.pngi
  cellcount.list[[paste0("Res.",res)]] <- table(VIC.data$Sample,VIC.data@meta.data[[paste0("integrated_snn_res.",res)]])
}
write.xlsx(cellcount.list,file="VIC_ClusterCount.xlsx",rowNames =TRUE)

#VIC老师确定分辨率选择0.2
Idents(VIC.data) <- "integrated_snn_res.0.2"
VIC.data <- RenameIdents(VIC.data,
                         '0'="VIC0",
                         '1'="VIC1",
                         '2'="VIC2",
                         '3'="VIC3",
                         '4'="VIC4",
                         '5'="VIC5")
VIC.data$sub.anno <- VIC.data@active.ident
saveRDS(VIC.data,"./../VIC/VIC.Cells_Integration_UMAP_annotation.rds")
#####################################################################
#GO-BP and KEGG analysis
setwd('./Analysis/week20220704/data/VIC')
vic.deg.raw <- read.csv("VIC_res.0.2_Find_AllMarkers.tab",sep="\t")
head(vic.deg.raw,3)

vic.deg.raw %>% filter(avg_logFC >0.5) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> VIC.top50
VIC.data <- ScaleData(VIC.data)
pdf('VIC_subcluster_heatmapDEGs.pdf')
DoHeatmap(VIC.data, features = VIC.top50$gene,group.by = "sub.anno",label=FALSE)+theme(axis.text.y = element_text(size = 2))
dev.off()
##########################################
#enrichment analysis
vic.deg.raw %>% filter(pct.1 >0.25) %>% filter(p_val < 0.05) -> vic.deg.filtered
anno <- c("VIC0","VIC1","VIC2","VIC3","VIC4","VIC5")
names(anno) <- c("0","1","2","3","4","5")
vic.deg.filtered$anno <- anno[as.character(vic.deg.filtered$cluster)]
VIC.enrich <- enrich_fuwai_new(degs=vic.deg.filtered,OrgDb="org.Hs.eg.db",organism="hsa",prefix="VIC")
p2 <- clusterProfiler::dotplot(VIC.enrich$kegg$`logfc>0.5`, showCategory=5, includeAll=FALSE,label_format=100)+
  labs(x='',y='',title="KEGG Enrichment Analysis")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",colour = "black",size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        plot.margin = margin(t = 20, b = 20, r = 30,  l = 30))
ggsave("VIC_KEGG_logfc>0.5_top5Dotplot.pdf",p2,width = 8,height = 8)

p2 <- clusterProfiler::dotplot(VIC.enrich$go$`logfc>0.75`, showCategory=5, includeAll=FALSE,label_format=100)+
  labs(x='',y='',title="GO Enrichment Analysis")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",colour = "black",size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        plot.margin = margin(t = 20, b = 20, r = 20,  l = 20))
ggsave("VIC_GOBP_logfc>0.75_top5Dotplot.pdf",p2,width = 9,height = 9)

p2 <- clusterProfiler::dotplot(VIC.enrich$kegg$`logfc>0.75`, showCategory=5, includeAll=FALSE,label_format=100)+
  labs(x='',y='',title="KEGG Enrichment Analysis")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",colour = "black",size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        plot.margin = margin(t = 10, b = 10, r = 80,  l = 80))
ggsave("VIC_KEGG_logfc>0.75_top5Dotplot.pdf",p2,width = 9,height = 9)

p2 <- clusterProfiler::dotplot(VIC.enrich$go$`logfc>1`, showCategory=5, includeAll=FALSE,label_format=100)+
  labs(x='',y='',title="GO Enrichment Analysis")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",colour = "black",size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        plot.margin = margin(t = 20, b = 20, r = 10,  l = 10))
ggsave("VIC_GOBP_logfc>1_top5Dotplot.pdf",p2,width = 9,height = 9)

p2 <- clusterProfiler::dotplot(VIC.enrich$kegg$`logfc>1`, showCategory=5, includeAll=FALSE,label_format=100)+
  labs(x='',y='',title="KEGG Enrichment Analysis")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",colour = "black",size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        plot.margin = margin(t = 10, b = 10, r = 100,  l = 100))
ggsave("VIC_KEGG_logfc>1_top5Dotplot.pdf",p2,width = 9,height = 9)

######################################################################
#SCENIC
#输出pySCINIC分析的输入文件
VIC.exprMat <- t(as.matrix(VIC.data@assays$RNA@counts))
write.table(VIC.exprMat,"VIC_RawCounts.csv",sep="\t",quote=F)

#pySCENIC运行完后绘图
setwd("./Analysis/week20220704/data/VIC/scenic")
write.table(data.frame(VIC.data$sub.anno),"VIC.celltype",sep="\t",quote=F)
vic.auc.cluster <- openxlsx::read.xlsx("VIC_auc_cluster.xlsx",rowNames = TRUE,colNames=TRUE)
vic.auc.scale <- as.data.frame(t(scale(vic.auc.cluster,center=TRUE,scale=TRUE)))
vic.auc.cluster.t <- data.frame(t(vic.auc.cluster))

diff.auc <- gsva.diff(vic.auc.cluster.t)
diff.auc %>% group_by(cluster) %>% top_n(n = 20, wt = logFC) -> vic.top20.regulons
#vic.top20.regulons <- data.frame(apply(vic.auc.cluster.t, 2, function(x){
#  rownames(head(vic.auc.cluster.t[order(x,decreasing = TRUE),],20))}))
diff.auc %>% group_by(cluster) %>% top_n(n = 3, wt = logFC) ->vic.top3.regulons 

vic.auc.cell <- read.csv("auc_mtx.csv",row.names = 1, check.names = FALSE)
VIC.data[["auc"]] <- NULL #删掉assay
VIC.data[["auc"]] <- CreateAssayObject(data=t(vic.auc.cell))

tf.regulons <- vic.auc.scale[unique(setorderv(vic.top3.regulons,c("cluster","logFC"),c(1,-1))$pathway),levels(VIC.data$sub.anno)]
#tf.regulons <- vic.auc.scale[unique(top3.tf),levels(VIC.data$sub.anno)]
#rownames(tf.regulons) <- regulons[rownames(tf.regulons),]
h.a <- pheatmap::pheatmap(head(vic.auc.scale,15), #fontsize_row=3, 
                   color=colorRampPalette(c("darkblue","white", high="darkred"))(100),
                   treeheight_row=10, treeheight_col=10, border_color=NA,
                   cellwidth = 20, cellheight = 12,
                   cluster_cols = FALSE, cluster_rows = FALSE,
                   legend_breaks = c(min(head(vic.auc.scale,15)), 0, max(head(vic.auc.scale,15))),
                   legend_labels = c("Low", "" ,"High"))

tf.express <- AverageExpression(VIC.data,group.by = "sub.anno",assays="RNA",
                                return.seurat = TRUE,
                                features = stringr::str_replace(unique(setorderv(vic.top3.regulons,c("cluster","logFC"),c(1,-1))$pathway), "\\(\\+\\)$", ""))

tf.data <- tf.express@assays$RNA@scale.data
h.tf <- pheatmap::pheatmap(tf.data, #fontsize_row=3, 
                   color=colorRampPalette(c("darkblue","white", high="darkred"))(100),
                   treeheight_row=10, treeheight_col=10, border_color=NA,
                   cellwidth = 20, cellheight = 12,
                   cluster_cols = FALSE,cluster_rows = FALSE,
                   legend_breaks = c(min(tf.data), 0, max(tf.data)),
                   legend_labels = c("Low", "", "High"))
pdf('test.pdf',h=6,w=6)
print(h.regulon)
print(h.a)
dev.off()
scenic.plot(adata=VIC.data,prefix="VIC")

#############################################################################
#cell cycle score
#Seurat CellCycleScoring
setwd("./Analysis/week20220704/data/VIC")
s.genes <- CaseMatch(search=cc.genes$s.genes ,match=rownames(VIC.data))
g2m.genes <- CaseMatch(search=cc.genes$g2m.genes,match=rownames(VIC.data))
set.seed(5)
VIC.cyclescore <- CellCycleScoring(VIC.data, s.features = s.genes, g2m.features = g2m.genes)
FeaturePlot(VIC.cyclescore,features = c("S.Score","G2M.Score"),cols=c("grey","Red1","Green"),blend = TRUE, blend.threshold = 0.00001,pt.size=0.5)
saveRDS(VIC.cyclescore,'VIC_Umap_anno_CellCycleScore.rds')
score.data <- data.frame(VIC.cyclescore$sub.anno,VIC.cyclescore$S.Score,VIC.cyclescore$G2M.Score,VIC.cyclescore$Phase)
colnames(score.data) <- c("Cluster","S.Score","G2M.Score","Phase")
write.xlsx(score.data,'VIC_CellCycleScores.xlsx',rowNames =TRUE)
#scran包的cyclone函数计算
#BiocManager::install("scran")
library("scran")
assignments <- readRDS('./cecllcycle/VIC_cyclone_CellCycle.rds')
anno.df <- VIC.data@meta.data[,c('Sample','Group','sub.anno')]
anno.df$Phases <- assignments$phases
vic.scores <- data.frame(t(assignments$scores))
colnames(vic.scores) <- rownames(anno.df)
pheatmap(vic.scores,
         show_colnames = F,
         #cluster_cols = FALSE, 
         cluster_rows = FALSE,
         annotation_col=anno.df,
         filename="VIC_CellScore_pheatmap.pdf",
         width=15,height = 8,
         cellwidth = 0.02, cellheight = 130)

################################################################################
#Group 差异分析
setwd('./Analysis/week20220718/data/VIC')
#Rcpp::sourceCpp("/BGFS1/home/kangbx/projects/decpp/fast_de.cpp")
#vic.markers.raw<-t_test_1vA(VIC.data@assays$RNA@data,as.character(VIC.data@meta.data[,"Group"]))
Idents(VIC.data) <- 'Group'
#vic.markers.raw <- FindAllMarkers(VIC.data,test.use = "t", min.cells.group = 0, min.pct=0, logfc.threshold = -Inf)
filter.markers <- function(markers=markers,prefix=prefix){
  markers %>% filter(pct.1 >0.3) %>% group_by(cluster) -> markers_flt_pct0.3
  markers_flt_pct0.3 %>% split(markers_flt_pct0.3$cluster) -> split.cluster
  for (i in 1:length(split.cluster)){
    library(data.table)
    split.cluster[[i]] <- setorderv(split.cluster[[i]],c("avg_logFC","pct.1","pct.2"),c(-1,1,-1))
  }
  openxlsx::write.xlsx(data.table::rbindlist(split.cluster),file=paste0(prefix,"_Find_AllMarkers_filtered.xlsx"), rowNames=FALSE)
  write.table(markers,file=paste0(prefix,"_Find_AllMarkers.tab"),quote=FALSE,sep="\t")
}
filter.markers(markers=vic.markers.raw,prefix="VIC.Group")

valcano.plot(vic.markers.raw,prefix="VIC")
vic.markers.raw %>% filter(pct.1 >0.3) %>% filter(avg_logFC>0) %>% 
  group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) -> vic.group.top100
enrich(markers = vic.group.top100,prefix='VIC.Group',OrgDb="org.Hs.eg.db",organism="hsa")
#######################################################################################
###GSEA分析需用FindAllMarkers计算差异基因
Idents(seurat.data) <- 'Group'
markers.raw <- FindAllMarkers(seurat.data,test.use = "t", min.cells.group = 0, min.pct=0, logfc.threshold = -Inf)
setwd('./Analysis/week20220718/data/VIC')
vic.markers.raw <- readRDS('VIC.Group_Allmarkers.rds')
fuwai.deg.enrich(markers.raw = vic.markers.raw, prefix = 'VIC')

#########################################################################################
#组间GSEA分析
#BiocManager::install("ReactomePA")
#install.packages("msigdbr")
library(ReactomePA)
library(msigdbr)
gsea.markers <- vic.markers.raw[vic.markers.raw$cluster=="MR",]
gsea.markers %>% filter(p_val < 0.05, avg_logFC > log(1.5)|avg_logFC < -log(1.5)) -> gsea.markers
gene.id<- bitr(gsea.markers$gene, fromType="SYMBOL",toType="ENTREZID",OrgDb=OrgDb)
gsea.markers <- gsea.markers[gsea.markers$gene %in% gene.id$SYMBOL,]
gsea.markers.merge <- merge(gsea.markers,gene.id,by.x='gene',by.y='SYMBOL')
geneList <- gsea.markers.merge$avg_logFC
names(geneList) <- as.character(gsea.markers.merge$ENTREZID)
geneList <- sort(geneList, decreasing = TRUE)

gsea.xlsx <- list()
#reactome
reactome.gsea <- ReactomePA::gsePathway(geneList, 
                                        organism = "human",
                                        pvalueCutoff = 1,
                                        pAdjustMethod = "BH", 
                                        verbose = FALSE,
                                        by='fgsea')
reactome.gsea <- setReadable(reactome.gsea,OrgDb,keyType="ENTREZID")
reactome.gsea@result <- reactome.gsea@result[reactome.gsea@result$pvalue <0.05,]
gsea.xlsx[['Reactome.GSEA']] <- data.frame(reactome.gsea)
data.gsea <- reactome.gsea
scores <- data.gsea@result$NES
names(scores) <- c(seq(1:length(scores)))
scores <- sort(scores,decreasing = TRUE)
MR.scores.index <- names(head(scores[scores>0],20))
NC.scores.index <- names(tail(scores[scores<0],20))
plot.list <- plot.gseaplot2(data.gsea = data.gsea,index.list=MR.scores.index)
cowplot::plot_grid(plotlist = plot.list,ncol=3)
ggsave(file=paste0("VIC.MR_Reactome_gseaplot2.pdf"),width = 12,height = 7)

data.gsea <- enrichplot::pairwise_termsim(data.gsea)
####hallmark GSEA
library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
hallmark.gsea <- GSEA(geneList, TERM2GENE=m_t2g, pvalueCutoff=1)
hallmark.gsea <- setReadable(hallmark.gsea,OrgDb,keyType="ENTREZID")
####GO GSEA
go.gsea <- gseGO(geneList = geneList,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              pvalueCutoff = 1)
go.gsea <- setReadable(hallmark.gsea,OrgDb,keyType="ENTREZID")
#################################################################
#GSEA result
gsea.result <- GSEA.Analysis(markers.raw=vic.markers.raw,group='MR',prefix="VIC",OrgDb="org.Hs.eg.db",organism="hsa")
#reactome plot
gsea.data <- gsea.result[['reactome']]
scores <- gsea.data@result$NES
names(scores) <- c(seq(1:length(scores)))
scores <- sort(scores,decreasing = TRUE)
MR.scores.index <- names(head(scores[scores>0],20))
NC.scores.index <- names(tail(scores[scores<0],20))
plot.list <- plot.gseaplot2(data.gsea = gsea.data,index.list=MR.scores.index)
cowplot::plot_grid(plotlist = plot.list,ncol=3)
ggsave(file=paste0("VIC.MR_Reactome_gseaplot2.pdf"),width = 12,height = 7)

gsea.data <- enrichplot::pairwise_termsim(gsea.data)
#hallmark plot
gsea.data <- gsea.result[['hallmark']]
scores <- gsea.data@result$NES
names(scores) <- c(seq(1:length(scores)))
scores <- sort(scores,decreasing = TRUE)
MR.scores.index <- names(head(scores[scores>0],20))
NC.scores.index <- names(tail(scores[scores<0],20))
plot.list <- plot.gseaplot2(data.gsea = gsea.data,index.list=MR.scores.index)
cowplot::plot_grid(plotlist = plot.list,ncol=2)
ggsave(file=paste0("VIC.MR_hallmark_gseaplot2.pdf"),width = 8,height = 3)

gsea.data <- enrichplot::pairwise_termsim(gsea.data)
#GO plot
gsea.data <- gsea.result[['go']]
scores <- gsea.data@result$NES
names(scores) <- c(seq(1:length(scores)))
scores <- sort(scores,decreasing = TRUE)
MR.scores.index <- names(head(scores[scores>0],20))
NC.scores.index <- names(tail(scores[scores<0],20))
plot.list <- plot.gseaplot2(data.gsea = gsea.data,index.list=MR.scores.index)
pdf(file=paste0("VIC.MR_GO_gseaplot2.pdf"),width = 20,height = 10)
print(cowplot::plot_grid(plotlist = plot.list[1:12],ncol=4))
print(cowplot::plot_grid(plotlist = plot.list[13:20],ncol=4,nrow=3))
dev.off()
plot.list <- plot.gseaplot2(data.gsea = gsea.data,index.list=NC.scores.index)
cowplot::plot_grid(plotlist = plot.list,ncol=2)
ggsave(file=paste0("VIC.NC_GO_gseaplot2.pdf"),width = 10,height = 8)

gsea.data <- enrichplot::pairwise_termsim(gsea.data)
#install.packages("ggnewscale")
#pdf('VIC_Reactome_emapplot.pdf')
##########################################################################################
#group SCENIC
setwd('./Analysis/week20220718/data/VIC/scenic')
write.table(data.frame(VIC.data$Group),"VIC.group.celltype",sep="\t",quote=F)
#/BGFS1/home/yanyl/miniconda3/envs/pyscenic/bin/python ./AUC_CellToCluster.py --in_dir . --out_dir . --prefix VIC.group
auc.data <- VIC.data
auc.cell <- read.csv("auc_mtx.csv",row.names = 1, check.names = FALSE)
if ("auc" %in% Assays(auc.data)){
  auc.data[["auc"]] <- NULL #删掉assay
}
auc.data[["auc"]] <- CreateAssayObject(data=t(auc.cell))
prefix <- "VIC.group"
#Regulon Specificity Score (RSS)
rss.cluster <- openxlsx::read.xlsx(paste0(prefix,"_rss.xlsx"),rowNames = TRUE,colNames=TRUE)
rss <- data.frame(t(rss.cluster))
top_regulons <- list()
rss %>% arrange(desc(MR)) %>% head(20) %>% rownames() -> top_regulons[['MR']]
rss %>% arrange(desc(NC)) %>% head(20) %>% rownames() -> top_regulons[['NC']]

target.genes <- read.csv('target_genes.xls',sep='\t')
gene.list <- c()
for (t.genes in target.genes$Target.Genes){
  for (gene in strsplit(t.genes,split=",")){
    gene.list <- c(gene.list,gene)
  }
}
gene.list <- gene.list[!duplicated(gene.list)]
DefaultAssay(auc.data) <- 'RNA'
Idents(auc.data) <- 'Group'
target.deg <- FindAllMarkers(auc.data, features = gene.list, test.use = "t", min.pct = 0.1, logfc.threshold = 0.25)
sort.markers <- function(markers=markers,prefix=prefix){
  markers %>% split(markers$cluster) -> split.cluster
  for (i in 1:length(split.cluster)){
    library(data.table)
    split.cluster[[i]] <- setorderv(split.cluster[[i]],c("avg_log2FC","pct.1","pct.2"),c(-1,1,-1))
  }
  openxlsx::write.xlsx(data.table::rbindlist(split.cluster),file=paste0(prefix,"_TargetGenes_DEG_filtered.xlsx"), rowNames=FALSE)
}
sort.markers(markers=target.deg,prefix=prefix)
auc.data <- ScaleData(auc.data,features = unique(target.deg$gene))
setorderv(target.deg,c("avg_log2FC"),c(-1))
heat.gene <-target.deg[target.deg$cluster=='MR',]$gene
pdf('VIC_TargetGenes_DEG_heatmap.pdf')
DoHeatmap(auc.data,
          features = heat.gene, 
          group.by = 'Group',
          label=FALSE)+
  theme(axis.text.y = element_text(size = 2))
dev.off()
saveRDS(auc.data,'VIC_addSCENIC_GroupRSS.rds')

tf.express <- AverageExpression(auc.data,group.by = "Group", assays="RNA",features = target.genes$TF )
target.deg.express <- AverageExpression(auc.data,group.by = "Group", assays="RNA",features = unique(target.deg$gene))
write.xlsx(list("TF"=data.frame(tf.express$RNA),"Target.Genes"=target.deg.express$RNA),file="VIC.group_TF&Target.Genes_AverageExpression.xlsx",rowNames=T)
if (F){
saveRDS(unique(target.deg$gene),'target.deg.rds')
anno.group <- subset(auc.data@meta.data, select = "Group")

pheatmap::pheatmap(auc.data@assays$RNA@scale.data[target.deg[target.deg$cluster=='MR',]$gene,], #fontsize_row=3, 
                   color=colorRampPalette(c("darkblue","white", high="darkred"))(100),
                   cellwidth = 0.1, cellheight = 1,
                   treeheight_row=10, treeheight_col=10, 
                   border_color=NA,
                   annotation_col = anno.group,
                   cluster_cols = FALSE,
                   #cluster_rows = FALSE,
                   legend_breaks = c(min(auc.data@assays$RNA@scale.data), 0, max(auc.data@assays$RNA@scale.data)),
                   legend_labels = c("Low", "", "High"),
                   main="Target_DEGs",
                   filename=paste0(prefix,"_Target_DEGs_heatmap.pdf"),
                   width=6,height = 6)

pheatmap::pheatmap(target.deg.express$RNA, #fontsize_row=3, 
                   color=colorRampPalette(c("darkblue","white", high="darkred"))(100),
                   cellwidth = 60, cellheight = 1,
                   treeheight_row=10, treeheight_col=10, 
                   border_color=NA,
                   #annotation_col = anno.group,
                   #cluster_cols = FALSE,cluster_rows = FALSE,
                   #legend_breaks = c(min(target.deg.express$RNA), 0, max(target.deg.express$RNA)),
                   legend_labels = c("Low", "", "High"),
                   main="Target_DEGs",
                   filename=paste0(prefix,"test.pdf"),
                   width=6,height = 10)
}
##################################################################