
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
library(cowplot)

setwd(".//Analysis/week20220613/data")
OutputDir = "./"
InputDir = "./"
Species = "hs"
FullLoad=TRUE
MultiCR=FALSE
AutoBasic<-function(OutputDir = "./",
                    InputDir = "./",
                    Species = "hs",
                    FullLoad = TRUE,
                    MultiCR = TRUE){
  
  pb <- txtProgressBar(style=3)
  star_time <- Sys.time()
  time <- 19
  
  ################################################################################
  #Parameter Configuration########################################################
  if (Species == "mm") {
    if (FullLoad) {
      OrgDb <- org.Mm.eg.db
    }
    MitochondrialPattern <- "^mt-"
    Angrycell_Species <- "mm"
  }else{
    if (FullLoad) {
      OrgDb <- org.Hs.eg.db
    }
    MitochondrialPattern <- "^MT-"
    Angrycell_Species = "hs"
  }
  print("2. Parameter Configuration")
  setTxtProgressBar(pb, 2/time)
  ################################################################################
  #Data Locked####################################################################
  DirList <- list.files(InputDir)
  #-------------------------------------------------------------------------------
  if (MultiCR) {
    GEXList <- c()
    for (i in DirList) {
      if (file.exists(paste0(InputDir,"/",i,"/outs/count"))) {
        GEXList <- c(GEXList,i)
      }
    }
    print(paste0("Total ",length(GEXList)," Transcriptome file(s) Locked"))
  }else{
    GEXList <- c()
    DirList <- list.files(paste0(InputDir,"/1_GEX"))
    for (i in DirList) {
      if (file.exists(paste0(InputDir,"/1_GEX/",i,"/outs"))){
        GEXList <- c(GEXList,i)
      }
    }
    print(paste0("Total ",length(GEXList)," Transcriptome file(s) Locked"))
  }
  ##Directory Create#-------------------------------------------------------------
  dir.create(paste0(OutputDir,"/Rdata"))
  dir.create(paste0(OutputDir,"/Stat"))
  dir.create(paste0(OutputDir,"/Img"))
  dir.create(paste0(OutputDir,"/Result"))
  
  if (FullLoad) {
    dir.create(paste0(OutputDir,"/Process"))
  }
  
  rm(DirList)
  print("3. Data Locked")
  setTxtProgressBar(pb, 3/time)
  ################################################################################
  #Sequence quality###############################################################
  if (length(GEXList)>0) {
    MetricsSummaryBook <- createWorkbook()
    for (i in GEXList) {
      if(MultiCR){
        MetricsSummaryTmp <- read_csv(paste0(InputDir,"/",i,"/outs/count/metrics_summary.csv"))
      }else{
        MetricsSummaryTmp <- read_csv(paste0(InputDir,"/1_GEX/",i,"/outs/metrics_summary.csv"))
      }
      addWorksheet(MetricsSummaryBook, paste0("Sample_",i))
      writeData(MetricsSummaryBook, paste0("Sample_",i),MetricsSummaryTmp)
    }
    saveWorkbook(MetricsSummaryBook, paste0(OutputDir,"/Stat/MetricsSummaryBook.xlsx"))
  }
  print("4. Sequence Quality Summary Done")
  setTxtProgressBar(pb, 4/time)
  rm(MetricsSummaryTmp,MetricsSummaryBook)
  ################################################################################
  #Load Data######################################################################
  for (j in GEXList) {
    i <- paste0("Data_",j)
    if(MultiCR){
      suerat.data<- Read10X(data.dir = paste0(InputDir,"/",j,"/outs/count/filtered_feature_bc_matrix"))
    }else{
      suerat.data<- Read10X(data.dir = paste0(InputDir,"/1_GEX/",j,"/outs/filtered_feature_bc_matrix"))
    }
    suerat_data<- CreateSeuratObject(counts = suerat.data, 
                                     project = i)
    assign(i,suerat_data)
  }
  if (length(GEXList)>1) {
    other <- ls(pattern="Data_")[-1]
    m <- function(x){eval(as.name(x))}
    data <- merge(x = get(ls(pattern="Data_")[1]),
                  y = c(lapply(X = other, function(x) m(x))),
                  add.cell.ids = eval(ls(pattern="Data_")),
                  project = "Data")
  }else{
    data <- suerat_data
  }
  saveRDS(data,file = paste0(OutputDir,"/Rdata/0-Data_RawObject.rds"))
  print(paste0("5. Load Data, Total ",length(GEXList)," file(s)"))
  setTxtProgressBar(pb, 6/time)
  rm(suerat_data,suerat.data,other,m,i,j)
  ################################################################################
  #Load Data######################################################################
  data[["percent.mt"]] <- PercentageFeatureSet(object = data, pattern = MitochondrialPattern)
  #data$Sample <- gsub("Data_","",data$orig.ident)
  Idents(data) <- "orig.ident"
  data <- RenameIdents(data,
                       "Data_201579Q_V24_M"="NC1",
                       "Data_201579AG_V49_Mit"="NC2",
                       "Data_201579AN_V59_Mit"="NC3",
                       "Data_201579E_V3_Mit"="MR1",
                       "Data_201579AA_V38_M"="MR2",
                       "Data_201579AD_V40_M"="MR3")
  data$Sample <- data@active.ident
  ###################
  #筛选表达HBB的细胞
  HBB.data <- data@assays$RNA@data["HBB",]
  HBB.Cells <- names(HBB.data[HBB.data > 0])
  length(HBB.Cells)
  data$HBB.state <- ifelse(rownames(data@meta.data) %in% HBB.Cells, "HBB.Pos","HBB.Neg")
  nc.sample <- c("NC1","NC2","NC3")
  data$Group <- ifelse(data$Sample %in% nc.sample,"NC","MR")
  saveRDS(data,file = paste0(OutputDir,"/Rdata/0-Data_RawObject.rds"))
  raw.data <- data
  head(HBB.Cells)
  CellCount_HBB <- as.data.frame(table(data$HBB.state,data$Sample)["HBB.Pos",])
  HBB.CellPercentage.Raw <- as.data.frame(prop.table(table(data$HBB.state,data$Sample))["HBB.Pos",])
  CellCount_Raw <- as.data.frame(table(data$Sample))
  UMICountMedian_Raw <- as.data.frame(tapply(data$nCount_RNA, INDEX=data$Sample, FUN=median))
  UMICountMean_Raw <- as.data.frame(tapply(data$nCount_RNA, INDEX=data$Sample, FUN=mean))
  FeatureMedian_Raw <- as.data.frame(tapply(data$nFeature_RNA, INDEX=data$Sample, FUN=median))
  FeatureMean_Raw <- as.data.frame(tapply(data$nFeature_RNA, INDEX=data$Sample, FUN=mean))
  MitochondrialPercentageMedian_Raw <- as.data.frame(tapply(data$percent.mt, INDEX=data$Sample, FUN=median))
  MitochondrialPercentageMean_Raw <- as.data.frame(tapply(data$percent.mt, INDEX=data$Sample, FUN=mean))
  P_Vln1 <- VlnPlot(object = data.raw, 
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    ncol = 3, 
                    pt.size = 0)
  #-------------------------------------------------------------------------------
  
  data <- subset(x = data, subset = percent.mt < 10 & nFeature_RNA > 200 & nFeature_RNA <4000)
  data <- subset(x=data,HBB.state=="HBB.Neg")
  
  CellCount <- as.data.frame(table(data$Sample))
  UMICountMedian <- as.data.frame(tapply(data$nCount_RNA, INDEX=data$Sample, FUN=median))
  UMICountMean <- as.data.frame(tapply(data$nCount_RNA, INDEX=data$Sample, FUN=mean))
  FeatureMedian <- as.data.frame(tapply(data$nFeature_RNA, INDEX=data$Sample, FUN=median))
  FeatureMean <- as.data.frame(tapply(data$nFeature_RNA, INDEX=data$Sample, FUN=mean))
  MitochondrialPercentageMedian <- as.data.frame(tapply(data$percent.mt, INDEX=data$Sample, FUN=median))
  MitochondrialPercentageMean <- as.data.frame(tapply(data$percent.mt, INDEX=data$Sample, FUN=mean))
  P_Vln2 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  #-------------------------------------------------------------------------------
  SeqStat<-cbind(CellCount_Raw, 
                 CellCount_HBB, 
                 CellCount, 
                 UMICountMedian_Raw, UMICountMedian,
                 UMICountMean_Raw,UMICountMean,
                 FeatureMedian_Raw, FeatureMedian,
                 FeatureMean_Raw,FeatureMean,
                 MitochondrialPercentageMedian_Raw, 
                 MitochondrialPercentageMedian,
                 MitochondrialPercentageMean_Raw,
                 MitochondrialPercentageMean)
  colnames(SeqStat) <- c("Sample",
                         "CellCount.Raw","HBBPositiveCellCount.Raw","Noneed","CellCount.Filtered",
                         "MedianUMICount.Raw","MedianUMICount.Filtered",
                         "MeanUMICount.Raw","MeanUMICount.Filtered",
                         "MedianGeneCount.Raw","MedianGeneCount.Filtered",
                         "MeanGeneCount.Raw","MeanGeneCount.Filtered",
                         "MedianPercentMT.Raw","MedianPercentMT.filtered",
                         "MeanPercentMT.Raw","MeanPercentMT.filtered")
  DF_SeqStat <- subset(SeqStat, select = -c(Noneed))
  openxlsx::write.xlsx(DF_SeqStat,file = paste0(OutputDir,"/Stat/","QualityControlState.xlsx"))
  #-------------------------------------------------------------------------------
  pdf(paste0(OutputDir,"/Img/","QualityControlState.pdf"),width = 16,height = 9)
  print(P_Vln1)
  print(P_Vln2)
  dev.off()
  saveRDS(data,file = paste0(OutputDir,"/Rdata/","1-Data_QualityControl.rds"))
  print(paste0("6. QC Finished ",ncol(data)," Cells Survived"))
  setTxtProgressBar(pb, 7/time)
  rm(CellCount_Raw,CellCount, UMICountMedian_Raw, UMICountMedian, FeatureMedian_Raw, FeatureMedian, MitochondrialPercentageMedian_Raw, MitochondrialPercentageMedian,P_Vln1,P_Vln2,SeqStat)
  ################################################################################
  #Pre-processing#################################################################
  data <- NormalizeData(object = data, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000,
                        verbose = FALSE)
  data <- FindVariableFeatures(object = data, 
                               selection.method = 'vst',
                               nfeatures = 2000,
                               verbose = FALSE)
  data <- ScaleData(object = data,
                    vars.to.regress = "percent.mt",
                    verbose = FALSE)
  saveRDS(data,file = paste0(OutputDir,"/Rdata/","2-Data_Pre-processing.rds")) 
  #-------------------------------------------------------------------------------
  top20 <- head(VariableFeatures(data), 20)
  P_TmpPlot1 <- VariableFeaturePlot(data)
  P_HighVariableFeatures <- LabelPoints(plot = P_TmpPlot1, 
                                        points = top20, 
                                        repel = TRUE)
  pdf(paste0(OutputDir,"/Img/","HighVariableFeatures.pdf"),width = 16,height = 9)
  print(P_HighVariableFeatures)
  dev.off()
  #-------------------------------------------------------------------------------
  if (FullLoad) {
    HighVariableGene <- subset(as.data.frame(data@assays$RNA@meta.features),vst.variable == "TRUE")
    write.xlsx(HighVariableGene,file = paste0(OutputDir,"/Stat/","HighVariableGene.xlsx"),rowNames = TRUE)
    rm(P_TmpPlot1,HighVariableGene,top20)
  }else{
    rm(P_TmpPlot1,top20)
  }
  print("9. Pre-processing Finished")
  setTxtProgressBar(pb, 10/time)
  ################################################################################
  #PCA & nPCs Selected############################################################
  data <- RunPCA(object = data, features = VariableFeatures(object = data),verbose = F)
  P_PCA <- DimPlot(object = data, 
                   reduction = "pca",
                   group.by = "Sample") + coord_fixed(1:1)
  P_Elbow <- ElbowPlot(data)
  pdf(paste0(OutputDir,"/Img/","DimensionalReduction_PCA.pdf"),width = 16,height = 9)
  print(P_PCA)
  print(P_Elbow)
  dev.off()
  #-------------------------------------------------------------------------------
  x <- data
  pct<-x[["pca"]]@stdev/sum(x[["pca"]]@stdev)*100
  cumu<-cumsum(pct)
  co1<-which(cumu >80 & pct<5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  nPCs <- min(co1,co2)
  print(paste0("10. PCA & nPCs Selected Finished. ",nPCs," PCs Selected"))
  #setTxtProgressBar(pb, 11/time)
  rm(x,co1,co2,cumu,pct) 
  ###nPCs=21
  ################################################################################
  #Clustering & tSNE/UMAP#########################################################
  #nPCs=13
  data <- FindNeighbors(object = data, 
                        reduction = "pca", 
                        dims = 1:nPCs, 
                        verbose = F)
  data <- FindClusters(object = data, 
                       verbose = F,
                       resolution = c(0.3,0.5,0.7))
                       #resolution = c(seq(.1,1.2,.2)))
  p_clust <- clustree(data@meta.data, prefix = "RNA_snn_res.")
  pdf(file=paste0(OutputDir,"/Img/","Allcells_clustree_plot.pdf"),width = 16,height = 9)
  print(p_clust)
  dev.off()
  #data <- RunTSNE(object = data, 
  #                dims = 1:nPCs)
  data <- RunUMAP(object = data,
                  #umap.method=umap-learn,
                  dims = 1:nPCs)
  DF_ClusterStat <- as.data.frame.matrix(table(Idents(data), data$Sample))
  write.xlsx(DF_ClusterStat,file = paste0(OutputDir,"/Stat/","CellCountDistribution.xlsx"),rowNames = TRUE)
  saveRDS(data,file = paste0(OutputDir,"/Rdata/","3-Data_Umap.rds"))
}
p_cluster <-DimPlot(object = data,
                    reduction = "umap",
                    pt.size = 0.00001,
                    #label = TRUE,
                    label = FALSE,
                    raster=FALSE,
                    group.by="RNA_snn_res.0.5",
                    #split.by = "Sample",
                    label.size = 10,
                    #cols=colors.sample
                    #cols=colors.cluster
) + ggtitle(label = 'Umap') + coord_fixed(1:1)

P_sample <- DimPlot(object = data,
        reduction = "umap",
        pt.size = 0.00001,
        #label = TRUE,
        label = FALSE,
        raster=FALSE,
        group.by="Sample",
        #split.by = "Sample",
        label.size = 10,
        cols=colors.sample
        #cols=colors.cluster
) + ggtitle(label = 'Umap') + coord_fixed(1:1)

sample.split <- DimPlot(object = data,
                        reduction = "umap",
                        pt.size = 0.00001,
                        #label = TRUE,
                        label = FALSE,
                        raster=FALSE,
                        group.by="RNA_snn_res.0.5",
                        split.by = "Sample",
                        label.size = 10,
                        #cols=colors.sample
                        #cols=colors.cluster
) + ggtitle(label = 'Umap') + coord_fixed(1:1)
pdf('cluster_sample_plot.pdf',width = 16,height = 9)
print(p_cluster)
print(P_sample)
print(sample.split)
dev.off()

#######################################################################
#绘制过滤前质控图
adata <- data
colors.sample <- brewer.pal(9,"Set3")[1:6]
colors.group <- brewer.pal(8,"Accent")[1:2]
#nc.sample <- c("NC1","NC2","NC3")
#data.raw$Group <- ifelse(data.raw$Sample %in% nc.sample,"NC","MR")
VlnPlot(object = data, features ="HBB",ncol=1,pt.size = 1,cols=colors.sample)
QC_plot_filtered(data=adata,prefix = "FilteredData",colors.sample=colors.sample,colors.group=colors.group)
#######################################################################
#去除批次效应
#filtered.data <- subset(x = raw.data, subset = percent.mt < 10 & nFeature_RNA > 200 & nFeature_RNA <4000)
#filtered.data <- subset(x=filtered.data,HBB.state=="HBB.Neg")
data.list <- SplitObject(data,split.by = "Sample")
# normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list)
group.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)
# this command creates an 'integrated' data assay
data.combined <- IntegrateData(anchorset = group.anchors)
#Now we can run a single integrated analysis on all cells!
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(data.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, vars.to.regress = "percent.mt",verbose = FALSE)
data.combined <- RunPCA(data.combined, verbose = FALSE)
ElbowPlot(data.combined)#nPCs=18
#########################
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:nPCs)
data.combined <- FindClusters(data.combined, resolution = c(seq(.1,1.1,.2)))
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:nPCs)
DF_ClusterStat <- as.data.frame.matrix(table(data.combined$integrated_snn_res.0.5, data.combined$Sample))
write.xlsx(DF_ClusterStat,file = paste0(OutputDir,"/Stat/","CellCountDistribution.xlsx"),rowNames = TRUE)
saveRDS(data.combined,"./Rdata/4-Data_DeBatch-integrated.Umap.rds")
colors.cluster <- DiscretePalette(13)
QC_plot_filtered(data=data.combined,res="integrated_snn_res.0.5",prefix = "FilteredData1",colors.sample=colors.sample,colors.group=colors.group,colors.cluster = colors.cluster)
cluster1 <- DimPlot(object = data.combined,
        reduction = "umap",
        pt.size = 0.00001,
        label = TRUE,
        #label = FALSE,
        raster=FALSE,
        group.by="integrated_snn_res.0.5",
        #cols=colors.cluster
        #split.by = "Group",
        #ncol=7
        #label.size = 10,
        #cols=colors.sample
) + ggtitle(label = 'Umap') + coord_fixed(1:1)
#clustree(data@meta.data, prefix = "integrated_snn_res.")
cluster.group <- DimPlot(object = data.combined,
                         reduction = "umap",
                         pt.size = 0.00001,
                         #label = TRUE,
                         label = FALSE,
                         raster=FALSE,
                         group.by="Group",
                         #cols=colors.cluster
                         #split.by = "Group",
                         #ncol=7
                         #label.size = 10,
                         cols=colors.group
) + ggtitle(label = 'Umap_Group') + coord_fixed(1:1)

cluster.sample <- DimPlot(object = data.combined,
                         reduction = "umap",
                         pt.size = 0.00001,
                         #label = TRUE,
                         label = FALSE,
                         raster=FALSE,
                         group.by="Sample",
                         #cols=colors.cluster
                         #split.by = "Sample",
                         #ncol=7
                         #label.size = 10,
                         cols=colors.sample
) + ggtitle(label = 'Umap_Sample') + coord_fixed(1:1)
cluster.sample1 <- DimPlot(object = data.combined,
                          reduction = "umap",
                          pt.size = 0.00001,
                          #label = TRUE,
                          label = FALSE,
                          raster=FALSE,
                          group.by="integrated_snn_res.0.5",
                          #cols=colors.cluster
                          split.by = "Sample",
                          #ncol=7
                          #label.size = 10,
                          #cols=colors.sample
) + ggtitle(label = 'Umap_Sample') + coord_fixed(1:1)
pdf('Clusters_plot.pdf',h=8,w=14)
cluster.plot<-plot_grid(cluster1,cluster.group,cluster.sample,ncol=3)
print(cluster.plot)
print(cluster.sample1)
dev.off()


adata <- readRDS("./Rdata/3-Data_Umap.rds")

#####################################
#差异基因分析及可视化需要用"RNA"的assay
DefaultAssay(data.combined) <- "RNA"
marker.new <- read.xlsx(".//Analysis/week20220613/data/heart_markers.xlsx",sheet=1)
marker.new
gene.new <- c()
for (genes.line in marker.new$Markers){
  for (gene in strsplit(genes.line,split=";")){
    gene.new <- c(gene.new,gene)
  }
}
pdf('FeaturePlot_markers.pdf',h=15,w=9)
FeaturePlot(data.combined,features=gene.new,ncol=3,order=TRUE)
dev.off()
pdf('DotPlot_markers.pdf',h=6,w=12)
DotPlot2(data.combined,features=gene.new,group.by="integrated_snn_res.0.5", x.angle = 45, x.hjust = 0.5)+
  scale_colour_gradient2(low = "blue",mid = "white",high = "red")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
#DotPlot(data.combined,features=gene.new,group.by="integrated_snn_res.0.5",split.by = "Group",cols=c("blue","red","white"))
dev.off()
data.combined <- ScaleData(data.combined, features = gene.new, vars.to.regress = "percent.mt",verbose = FALSE)
pdf('test.pdf')
DoHeatmap(data.combined,features = gene.new,group.by="integrated_snn_res.0.5")
dev.off()
heart.top50 <- find.markers(data=data.combined,res="integrated_snn_res.0.5",cell.name="fuwai",n=50)
VlnPlot(data.combined,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "integrated_snn_res.0.7",pt.size=0)
VlnPlot(data.combined,features = gene.new,group.by = "integrated_snn_res.0.5",pt.size=0)
Idents(data.combined) <- "integrated_snn_res.0.5"
data.combined.anno <- RenameIdents(data.combined,
                                   "0"="VIC",
                                   "1"="VIC",
                                   "2"="VIC",
                                   "3"="LowQuality",
                                   "4"="VIC",
                                   "5"="VIC",
                                   "10"="VIC",
                                   "11"="VIC",
                                   "8"="VEC",
                                   "7"="Lymphocyte",
                                   "6"="Myeloid",
                                   "9"="Myeloid",
                                   "12"="Mast cell")
data.combined.anno$Annotation <- data.combined.anno@active.ident
filter.lowquality <- subset(data.combined.anno, Annotation!="LowQuality")
p1 <- DimPlot(object = filter.lowquality,
        reduction = "umap",
        pt.size = 0.00001,
        label = TRUE,
        #label = FALSE,
        raster=FALSE,
        #group.by="integrated_snn_res.0.5",
        #cols=colors.cluster
        #split.by = "Group",
        #ncol=7
        #label.size = 10,
        #cols=colors.sample
) + ggtitle(label = 'Umap') + coord_fixed(1:1)
p2 <- DimPlot(object = filter.lowquality,
              reduction = "umap",
              pt.size = 0.00001,
              label = TRUE,
              #label = FALSE,
              raster=FALSE,
              #group.by="integrated_snn_res.0.5",
              #cols=colors.cluster
              split.by = "Group",
              #ncol=7
              #label.size = 10,
              cols=color.cluster
) + ggtitle(label = 'Umap') + coord_fixed(1:1)

pdf('Celltype_plot.pdf',h=8,w=14)
print(p1)
print(p2)
dev.off()
rm(data.combined)
saveRDS(data.combined.anno,"./Rdata/4-Data_DeBatch-integrated.Umap.annotation.rds")
marker.new <- gene.new[(!gene.new %in% c("MYH7","MYL2","ACTC1"))] 

heat.png <- DoHeatmap(filter.lowquality,features =marker.new ,label=FALSE)
#color.cluster <- unique(ggplot_build(p1)$data[[1]]$colour) #获取p1绘图时使用的颜色
#color.cluster <- c("#F8766D","#A3A500","#00BF7D","#00B0F6")
cell.cluster <- unique(ggplot_build(p1)$data[[2]]$label)
names(cell.cluster) <- unique(ggplot_build(p1)$data[[1]]$colour)
color.cluster <- names(sort(cell.cluster))
vln.png <- vioplot2(filter.lowquality,
         paths=marker.new,
         facet.toward = "col",
         cell.order=rev(levels(ggplot_build(p1)$data[[2]]$label)),
         color.use = rev(color.cluster),
         remove.axis.x = T,
         panel.spacing = 0)
pdf('Heatmap_Plot.pdf',h=6,w=5)
print(heat.png)
dev.off()
pdf('Vlnplot_markers.pdf',h=3,w=9)
print(vln.png)
dev.off()
#VlnPlot(filter.lowquality,features=marker.new[1:3],pt.size = 0)

filter.lowquality$Annotation <- factor(filter.lowquality$Annotation)
group.compare.data <- as.data.frame(prop.table(table(filter.lowquality$Sample,filter.lowquality$Annotation),margin = 1))
group.compare.data$Sample <- rownames(group.compare.data)
group.compare.data$Group <- ifelse(group.compare.data$Var1 %in% c("MR1","MR2","MR3"),"MR","NC")
group.compare.data$ratio <- round(group.compare.data$Freq*100,1)
names(colors.sample) <- levels(group.compare.data$Var1)
color.use <- colors.sample[as.vector(group.compare.data$Var1)]
#sub.data <- group.compare.data[group.compare.data$Var2=="VIC",]

#比较组间差异显著性
ggpubr::compare_means(ratio~Group,group.compare.data,
                      method = "wilcox.test",paired = FALSE,
                      group.by = "Var2")

Cellcluster.top50 <- find.markers(filter.lowquality,res="Annotation",cell.name="CellCluster",n=50)
Cellcluster.top50 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> top5.genes
#调整基因顺序
top5.genes$cluster <- factor(top5.genes$cluster,levels=c("VIC","VEC","Lymphocyte","Myeloid","Mast cell"))
top5.genes <- top5.genes[order(top5.genes$cluster), ]
pdf('CellType_DotPlot_markers.pdf',h=5,w=10)
DotPlot2(filter.lowquality,features=top5.genes$gene,group.by="Annotation", x.angle = 45, x.hjust = 0.5)+
  scale_colour_gradient2(low = "blue",mid = "white",high = "red")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
#DotPlot(data.combined,features=gene.new,group.by="integrated_snn_res.0.5",split.by = "Group",cols=c("blue","red","white"))
dev.off()
saveRDS(filter.lowquality,"./Rdata/Integration.Umap.annotation_discardLowQuality.rds")

useful.cells <- filter.lowquality
useful.cells <- ScaleData(useful.cells, features = rownames(useful.cells))
top50 <- find.markers(data=useful.cells,res="Annotation",cell.name="all.markers",n=50)
top50$cluster <- factor(top50$cluster,levels = c("VIC","VEC","Lymphocyte","Myeloid","Mast cell"))
top50 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> top5
pdf('CellType_heatmap_top50.pdf')
DoHeatmap(useful.cells,features = top50[order(top50$cluster),]$gene, group.by="Annotation",label=FALSE)+theme(axis.text.y = element_text(size = 2))
dev.off()
