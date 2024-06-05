BiocManager::install("slingshot")
BiocManager::install("tradeSeq")
library(slingshot)
library(tradeSeq)
library(Seurat)
library(ggplot2)
library(dplyr)

#################################################################
#slingshot VIC
setwd('./week20220815/data/VIC')
Idents(VIC.data) <- 'sub.anno'
vic.slingshot <- subset(VIC.data, idents = c("VIC1", "VIC2", "VIC3"))
vic.crv <- slingshot(Embeddings(vic.slingshot, "umap"), vic.slingshot$sub.anno, start.clus = 'VIC3')
colors <- (scales::hue_pal())(6)[c(2,3,4)]
names(colors) <- c("VIC1", "VIC2", "VIC3")
saveRDS(list(sling=vic.crv, count=vic.slingshot@assays$RNA@counts),'VIC_slingshot.rds')
slingshot.plot(seurat.data=vic.slingshot,prefix='VIC',slingshot.data = vic.crv, colors=colors,method='umap')
vic.fitgam <- readRDS('VIC_fitgam.rds')
setwd('./pca')
sling.pca <- slingshot(Embeddings(vic.slingshot, "pca"), vic.slingshot$sub.anno, start.clus = 'VIC3')
slingshot.plot(seurat.data=vic.slingshot,prefix='VIC.pca',slingshot.data = sling.pca, colors=colors,method='pca')
saveRDS(list(sling=sling.pca, count=vic.slingshot@assays$RNA@counts),'VIC_slingshot_pca.rds')

top5.vic <- c('CCDC80','PCOLCE2','PI16','GSN','MFAP5',
              'FMOD','PRELP','CLU','ASPN','COMP',
              'EFEMP1','ADH1B','CA3','APOE','GGT5')  
vic.pca.fitgam <- readRDS('VIC.pca_fitgam.rds')
plot.pseudotime.gene(fitgam=vic.pca.fitgam,gene.count=vic.slingshot@assays$RNA@counts,prefix='VIC.pca',genes=top5.vic)
#######################################################################
#slingshot VEC
setwd('./week20220815/data/VEC')
vec.crv <- slingshot(Embeddings(VEC.data, "umap"), VEC.data$sub.anno, start.clus = 'VEC2')
all(colnames(VEC.data@assays$RNA@counts)==rownames(VEC.data$Sample))
vec.batch <- factor(VEC.data$Sample)
vec.U <- model.matrix(~vec.batch)
saveRDS(list(sling=vec.crv, count=VEC.data@assays$RNA@counts),'VEC_slingshot.rds')
colors <- (scales::hue_pal())(6)[c(1,2,3)]
names(colors) <- c("VEC0", "VEC1", "VEC2")
slingshot.plot(seurat.data=VEC.data,prefix='VEC',slingshot.data = vec.crv, colors=colors,method='umap')
top5.vec <- c('PIK3R3','IGFBP3','FABP4','A2M','CLDN5',
              'PLA2G2A','LEPR','COLEC11','NSG1','POSTN',
              'COL1A2','LUM','DCN','COL1A1','APOE')
vec.pca.fitgam <- readRDS('VEC_fitgam.rds')
plot.pseudotime.gene(fitgam=vec.pca.fitgam,gene.count=VEC.data@assays$RNA@counts,prefix='VEC',genes=top5.vec)

