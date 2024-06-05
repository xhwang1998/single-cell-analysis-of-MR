setwd('./Analysis/week20220725/data')
library(CellChat)
library(ggplot2)
library(Seurat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(dplyr)
library(grDevices)

useful.cell <- readRDS('Integration.Umap.annotation_discardLowQuality.rds')
object.list <- list()
for (group in unique(useful.cell$Group)){
  group.cells <- subset(useful.cell,Group==group)
  cellchat.group <- createCellChat(object = group.cells@assays$RNA@data,
                                   meta = group.cells@meta.data[,c('orig.ident','Sample','Group','Annotation')], 
                                   group.by = "Annotation")
  cellchat.group@DB <- CellChatDB.human
  cellchat.group %>% 
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb() %>%
    filterCommunication(min.cells = 10) %>%
    computeCommunProbPathway() %>%
    aggregateNet() -> object.list[[group]]
}
saveRDS(object.list,'Group_CellchatObjList_MR-NC.rds')
cellchat.obj <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat.obj,'Group_MergeCellchatObj.rds')
p1 <- compareInteractions(cellchat.obj, 
                    show.legend = F, 
                    group=c(1,2), 
                    measure = 'count',
                    xlabel = 'Differential number of interactions')+
  theme(plot.margin = margin(t=100,b=90,r=50,l=50),
        axis.title.x = element_text(size=11))
p2 <- compareInteractions(cellchat.obj, 
                          show.legend = F,
                          group=c(1,2), 
                          measure = "weight",
                          xlabel = 'Differential interaction strength')+
  theme(plot.margin = margin(t=100,b=90,r=50,l=50),
        axis.title.x = element_text(size=11))
# the second dataset compared to the first one
p3 <- netVisual_diffInteraction(cellchat.obj,
                                weight.scale = T,
                                measure = "count",
                                comparison = c(2,1),
                                margin = 0.1)
p4 <- netVisual_diffInteraction(cellchat.obj,
                                weight.scale = T,
                                measure = "weight",
                                comparison = c(2,1),
                                margin = 0.1)
pdf("CellChat_interaction_count&strength.pdf",width = 8,height = 6)
print(p1+p2)
print(p3)
print(p4)
dev.off()
p5 <- netVisual_heatmap(cellchat.obj,
                        measure = "count",
                        comparison = c(2, 1),
                        height = 4,
                        width=4)
p6 <- netVisual_heatmap(cellchat.obj,
                        measure = "weight",
                        comparison = c(2, 1),
                        height = 4,
                        width=4)

pdf("CellChat_interaction_count&strength_heatmap.pdf",width = 4.8,height = 4)
print(p5)
print(p6)
dev.off()
information.flow <- rankNet(cellchat.obj,
                            mode = "comparison",
                            stacked = T, 
                            do.stat = TRUE,
                            comparison = c(2, 1),
                            color.use = c("#00BFC4","#F8766D"),
                            return.data=TRUE)
pdf("CellChat_InformationFlow.pdf",width = 6,height = 10)
print(information.flow$gg.obj)
dev.off()