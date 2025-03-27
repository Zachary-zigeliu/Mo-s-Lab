library(vctrs)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(patchwork) 
library(igraph) 
library(dplyr)
library(future)
library(NMF)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(igraph)
library(dplyr)
library(cowplot)
library(openxlsx)
library(patchwork)
options(stringsAsFactors = FALSE)
dpi=300

########  CellChat  ########
Cell <- "CellChat"
dir = '~/Tissue/re1.4'  
folder <- paste0(dir,"/",Cell)
PURE_RDS <-"~/Tissue/re1.4/PURE_RDS"
setwd(dir);dir.create(Cell);setwd(folder)

##  -- Fibroblasts & T cells CellChat --
Fibro <- readRDS(paste0(PURE_RDS,"/Non_ImmuneCells/","Fibro","/PURE.RDS"))
Idents(Fibro) <- Fibro@meta.data[["new_seurat_clusters2"]]
Tcell <- readRDS(paste0(PURE_RDS,"/ImmuneCells/","Tcell","/PURE.RDS"))
Idents(Tcell) <- Tcell@meta.data[["new_seurat_clusters2"]]
Cells <- merge(Fibro,y= Tcell)
levels(Cells)

data.input <- GetAssayData(Cells, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Cells)
meta <- data.frame(labels = labels,row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels",levels = c(
  "Naive-IL7R",
  "CD4-FOXP3","CD4-CXCL13",
  "CD8-MKI67","CD8-STMN1","CD8-HAVCR2","CD8-HSPA1A","CD8-ISG15",
  "CD8-GNLY","CD8-GZMK",
  "F-HSD17B2","F-PTGDS","F-POSTN","F-SFRP2","F-PCOLCE2","F-IER3",
  "F-ALPL","F-CPE","F-TPPP3","F-ISG15","F-CXCL14","F-NGFR","F-STMN1")) # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

########  AllCellChatDB  通路数据库分析  ###########  
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat) 
cellchat <- aggregateNet(cellchat) 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot',type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot',type = "structural")
cellchat <- netClustering(cellchat, type = "structural")


########  CellChat  ########  作图 Figs
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
ggsave('scatter_AllDB.png', plot = gg1, width=5.5,height=5, dpi=300)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
png("outgoing_AllDB.png", width=dpi*7,height=dpi*15,units = "px",res = dpi,type='cairo')
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                         width = 12, height = 28,
                                         font.size = 10, font.size.title = 14)
ht1
dev.off()

#ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
png("incoming_AllDB.png", width=dpi*7,height=dpi*15,units = "px",res = dpi,type='cairo')
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                         width = 12, height = 28,
                                         font.size = 10, font.size.title = 14)
ht2
dev.off()


groupSize <- as.numeric(table(cellchat@idents))
levels(cellchat@idents) # 查看细胞类群数
png("Number_of_interactions.png", width=dpi*8,height=dpi*8,units = "px",res = dpi,type='cairo')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
png("Interaction_strength.png", width=dpi*8,height=dpi*8,units = "px",res = dpi,type='cairo')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

dir.create('pathways')
pathways <- cellchat@netP[["pathways"]]
for (pathways.show in pathways) {
  plot <- netAnalysis_contribution(cellchat, signaling = pathways.show)
  ggsave(paste0('pathways/',pathways.show,'.png'),plot = plot, width = 5, height = 5, dpi=300)
}

dir.create('pathways_network')
pathways <- cellchat@netP[["pathways"]]
for (pathways.show in pathways) {
  png(paste0('pathways_network/',pathways.show,'.png'), width=dpi*8,height=dpi*2.5,units = "px",res = dpi,type='cairo')
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 14, height = 2, font.size = 10)
  dev.off()
}

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
dir.create('Hierarchy plot')
vertex.receiver = c(1:10) # a numeric vector. 1-10是T细胞
# pathways.show <- c("CCL") 
pathways <- cellchat@netP[["pathways"]]
for (pathways.show in pathways) {
  png(paste0('Hierarchy plot/',pathways.show,'.png'), width=dpi*12.5,height=dpi*7.5,units = "px",res = dpi,type='cairo')
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
  dev.off()
}

## 使用小提琴/气泡图绘制信号基因表达分布
dir.create('VlnPlot')
lapply(pathways, function(i){
  gg1 <- plotGeneExpression(cellchat, signaling = i, angle.x = 45, colors.ggplot = T)
  ggsave(paste0("VlnPlot/",i,"_Vlin_plot.png"),plot = gg1, width = 7,height = 5, dpi = dpi)
})



