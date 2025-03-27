source('/data2/teacher/yang/R_package.R')

########  Fibro  ########
Cell <- "Fibro"
dir = "~/Tissue/re1.4"  ## 新服务器路径
folder <- paste0(dir,"/",Cell)
res <- "re0.9"

Integrated_1.4 <- readRDS("./Tissue/RDS/Tissue_Integrated_1.4.RDS")

Sub <- subset(Integrated_1.4, idents = c(0,5,14,25,33,37,42))  ## 提取总图中相应的细胞类群

##重新分类
Sub <- NormalizeData(Sub, normalization.method = "LogNormalize", scale.factor = 10000)
Sub <- FindVariableFeatures(Sub, selection.method = "vst", nfeatures = 2000)
Sub <- ScaleData(Sub, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE ) # "percent.mt","percent.ribo"
Sub <- RunPCA(Sub, npcs = 100, verbose = FALSE)

## Harmony去批次效应
Sub <- Sub %>% RunHarmony("orig.ident", plot_convergence = TRUE)

Sub <- RunUMAP(Sub, reduction = "harmony", dims = 1:30)
Sub <- FindNeighbors(Sub, reduction = "harmony", dims = 1:30)
Sub_0.9 <- FindClusters(Sub, resolution = 0.9)
Sub_0.9.MAST <- FindAllMarkers(Sub_0.9, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST")

## 去除混杂细胞后，提取相应的细胞类群
clus <- c(0,1,2,3,4,5,6,7,9,10,11,16,17)
PURE <- subset(Sub_0.9, idents = clus)  
PURE@meta.data[["seurat_clusters"]] <- droplevels(PURE@meta.data[["seurat_clusters"]])

## 命名 new_seurat_clusters
new.cluster.ids <- c('0','1','2','3','4','5','6','7','8','9','10','11','12')
names(new.cluster.ids) <- levels(PURE)
PURE <- RenameIdents(PURE, new.cluster.ids)
PURE@meta.data$new_seurat_clusters <- PURE@active.ident
PURE <- RunUMAP(PURE, reduction = "harmony", dims = 1:30)

## 命名 new_seurat_clusters2
new.cluster.ids <- c('F-PCOLCE2',
                     'F-SFRP2','F-PTGDS','F-CXCL14','F-POSTN','F-IER3',
                     'F-HSD17B2','F-CPE','F-ALPL','F-TPPP3','F-ISG15',
                     'F-NGFR','F-STMN1')
names(new.cluster.ids) <- levels(PURE)
PURE <- RenameIdents(PURE, new.cluster.ids)
PURE@meta.data$new_seurat_clusters2 <- PURE@active.ident

## 命名 new_seurat_clusters3
new.cluster.ids <-c('0    F-PCOLCE2',
                    '1    F-SFRP2',  '2    F-PTGDS','3    F-CXCL14','4    F-POSTN','5    F-IER3',
                    '6    F-HSD17B2','7    F-CPE',  '8    F-ALPL', '9    F-TPPP3','10  F-ISG15',
                    '11  F-NGFR',   '12  F-STMN1')
names(new.cluster.ids) <- levels(PURE)
PURE <- RenameIdents(PURE, new.cluster.ids)
PURE@meta.data$new_seurat_clusters3 <- PURE@active.ident

plot <- DimPlot(PURE, reduction = "umap", label = TRUE, label.size = 6)
ggsave("PURE_new3A.png", plot = plot, width = 7.5, height = 6, dpi=300)

Idents(PURE) <- PURE@meta.data[["new_seurat_clusters"]]
saveRDS(PURE,file = paste0(folder,"/RDS/PURE.RDS"))


library(RColorBrewer)
######  细胞比例
nb.cols <- length(levels(PURE@active.ident))
mycolors <- colorRampPalette(brewer.pal(nb.cols, "Set3"))(nb.cols)

cell.prop<-as.data.frame(prop.table(table(Idents(PURE), PURE$orig.ident_new4)))
colnames(cell.prop)<-c("Cluster","origin","Proportion")
cell.prop <- cell.prop[cell.prop$origin %in% c("Normal", "MIBC","NMIBC"),]  

plot <- ggplot(cell.prop,aes(origin,Proportion,fill=Cluster))+
  geom_bar(stat="identity",position="fill")+
  labs(x = '',y = 'Cell proportion', title = '')+
  theme(panel.background = element_rect(fill = NA,colour = 'black'),
        #panel.background = element_blank(),#"transparent"
        axis.ticks = element_blank(),#axis.line = element_line(color = 'black'),
        legend.text = element_text(size=14),
        axis.title =element_text(size = 20),
        axis.text.x =element_text(size = 14, color = 'black'),axis.text.y =element_text(size = 12, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 20))+
  scale_x_discrete(limits=c("Normal", "MIBC","NMIBC"))+   scale_y_continuous(labels = scales::percent)+ #纵坐标变为百分比
  scale_fill_manual(values = mycolors) + ##调用颜色方案
  RotatedAxis()#+
ggsave("CellProportion2.png", plot = plot, width = 4.5, height = 6, dpi=300)  


##  小提琴图   
PURE <- readRDS(paste0(folder,"/RDS/PURE.RDS"))
Idents(PURE) <- PURE@meta.data$new_seurat_clusters2
Idents(PURE) <- factor(Idents(PURE), 
                       levels =c('F-HSD17B2','F-PTGDS',
                                 'F-POSTN','F-SFRP2','F-PCOLCE2','F-IER3',
                                 'F-ALPL','F-CPE','F-TPPP3',
                                 'F-ISG15','F-CXCL14','F-NGFR','F-STMN1'))
markers <- c("HAS1", "HAS2")

for(marker in markers){
  plot <- VlnPlot(PURE, features = marker,pt.size = 0)+NoLegend()+
    xlab(element_blank()) + ylab(marker) + ggtitle(element_blank()) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),#legend.position = "none", 
          axis.title.y = element_text(size = rel(1), angle = 90, vjust = 0.5))
  ggsave(paste0("violin_",marker,".png"), plot = plot, width = 6, height = 1.25, dpi=300)
}


###   火山图
library(ggrepel)
dir.create("Volcano");setwd("Volcano")

Idents(PURE) <- PURE@meta.data[["orig.ident_new4"]] ## 设置 active.ident 为组织来源
tissue_markers <- FindMarkers(PURE, ident.1="MIBC", ident.2="NMIBC", logfc.threshold=0,verbose = TRUE) 
tissue_markers$gene <- rownames(tissue_markers)
PURE.logFC_cutoff <- with(tissue_markers,mean(abs(avg_log2FC)) + 2*sd(abs(avg_log2FC)))
tissue_markers$sig = as.factor(ifelse(tissue_markers$p_val_adj < 0.05 & abs(tissue_markers$avg_log2FC) > PURE.logFC_cutoff,ifelse(tissue_markers$avg_log2FC > PURE.logFC_cutoff ,'UP','DOWN'),'NO'))
tissue_markers$PURE.logFC_cutoff<-PURE.logFC_cutoff
gene_PURE <- c('MMP11','POSTN','COL12A1','COL5A2','MMP1','STAT1','ISG15','MMP13','MMP14',
               'MGP','IGF1','PLA2G2A','COL4A4','PCOLCE2')
PURE_selected <- tissue_markers[gene_PURE,]

p1<-ggplot(tissue_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +xlim(-1.5, 1.5) +
  geom_point(aes(color = sig), size=1) +labs(color="Significance",size=5)+#guides(fill=guide_legend(reverse=TRUE))+
  scale_color_manual(values = c("blue","grey", "red")) +
  geom_text_repel(data =  PURE_selected,aes(label = gene),  size = 4,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines")) +
  geom_vline(xintercept=c(-PURE.logFC_cutoff,PURE.logFC_cutoff),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)", y="-log10 (p-value)") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("volcano.png", plot = p1, width = 5, height = 5, dpi=300)


####         汽泡图
dir.create("DotPlot");setwd("DotPlot")
PURE <- readRDS(paste0(folder,"/RDS/PURE.RDS"))
Idents(PURE) <- PURE@meta.data$new_seurat_clusters2
Idents(PURE) <- factor(Idents(PURE), 
                       levels =c('F-PCOLCE2','F-CPE','F-IER3','F-ISG15',
                                 'F-CXCL14','F-SFRP2','F-ALPL','F-PTGDS','F-POSTN','F-HSD17B2','F-STMN1',
                                 'F-TPPP3','F-NGFR'))

markers <- c("TGFB2", "TGFB1", "TGFB3","IGF1","FGF7", "FGF2","FGF18","FGF9", #,"HGF","FGF10","IL33","IL6",
             "TGFBR1","TGFBR2","TGFBR3","FGFR1","IGF1R", "IGF2R") #, "FGFR2", "FGFR3", "FGFR4"
plot <- DotPlot(PURE, features = markers,cols = c("lightgrey","#FF0000"),col.min = 0,col.max = 2) + RotatedAxis() #c('white','#F8766D')
plot <- plot+scale_x_discrete(limits=markers) +labs(x=element_blank(), y=element_blank()) ### 调整cluster的顺序
plot[["guides"]][["size"]][["title"]] <- "Percent \nExpressed" ; plot[["guides"]][["colour"]][["title"]] <- "Average \nExpression"
plot
ggsave("Mark_GFs.png", plot = plot, width = 7, height = 4.75, dpi=300)

markers <- c("COL1A1", "COL1A2","COL6A1","COL6A2","COL6A3", "COL3A1","COL5A1","COL5A2",
             "COL12A1","COL16A1","COL8A1", "COL8A2","COL15A1","COL10A1", "COL11A1","COL7A1",
             "COL5A3","COL14A1","COL4A1","COL4A2", "COL4A3","COL4A4","COL24A1","COL4A5",
             "COL4A6","COL9A1","COL13A1")
Idents(PURE) <- factor(Idents(PURE), 
                       levels =c('F-NGFR','F-CPE','F-PCOLCE2','F-CXCL14','F-ISG15','F-IER3',
                                 'F-PTGDS','F-SFRP2','F-ALPL','F-POSTN','F-STMN1','F-TPPP3',
                                 'F-HSD17B2'))
plot <- DotPlot(PURE, features = markers,cols = c("lightgrey","#FF0000"),col.min = 0,col.max = 2) + RotatedAxis() #c('white','#F8766D')
plot <- plot+scale_x_discrete(limits=markers) +labs(x=element_blank(), y=element_blank()) ### 调整cluster的顺序
plot[["guides"]][["size"]][["title"]] <- "Percent \nExpressed" ; plot[["guides"]][["colour"]][["title"]] <- "Average \nExpression"
ggsave("Mark_COLs.png", plot = plot, width = 11, height = 4.75, dpi=300)

markers <- c("VEGFA","VEGFB","VEGFC","VEGFD" ,"PIGF","PROK1")
Idents(PURE) <- factor(Idents(PURE), 
                       levels =c('F-NGFR','F-TPPP3','F-CPE','F-PCOLCE2','F-ISG15','F-CXCL14','F-IER3',
                                 'F-PTGDS','F-SFRP2','F-ALPL','F-POSTN','F-STMN1',
                                 'F-HSD17B2'))
plot <- DotPlot(PURE, features = markers,cols = c("lightgrey","#FF0000"),col.min = 0,col.max = 2) + RotatedAxis() #c('white','#F8766D')
plot <- plot+scale_x_discrete(limits=markers) +labs(x=element_blank(), y=element_blank()) ### 调整cluster的顺序
plot[["guides"]][["size"]][["title"]] <- "Percent \nExpressed" ; plot[["guides"]][["colour"]][["title"]] <- "Average \nExpression"
ggsave("Mark_VEGFs.png", plot = plot, width = 5.25, height = 4.75, dpi=300)

markers <- c("MMP2","MMP14","MMP23B","MMP19")
Idents(PURE) <- factor(Idents(PURE), 
                       levels =c('F-NGFR','F-CPE','F-CXCL14','F-PCOLCE2','F-ISG15','F-IER3',
                                 'F-PTGDS','F-SFRP2','F-ALPL','F-POSTN','F-HSD17B2','F-STMN1','F-TPPP3'))
plot <- DotPlot(PURE, features = markers,cols = c("lightgrey","#FF0000"),col.min = 0,col.max = 2) + RotatedAxis() #c('white','#F8766D')
plot <- plot+scale_x_discrete(limits=markers) +labs(x=element_blank(), y=element_blank()) ### 调整cluster的顺序
plot[["guides"]][["size"]][["title"]] <- "Percent \nExpressed" ; plot[["guides"]][["colour"]][["title"]] <- "Average \nExpression"
ggsave("Mark_MMPs.png", plot = plot, width = 3.75, height = 4.75, dpi=300)

markers <- c("IFNA1","IFNA2","IFNB1","IFNG","IFNAR1","IFNAR2","IFNGR1","IFNGR2") #
Idents(PURE) <- factor(Idents(PURE), 
                       levels =c('F-STMN1','F-ALPL','F-PCOLCE2','F-TPPP3','F-CPE','F-ISG15','F-CXCL14','F-NGFR',
                                 'F-POSTN','F-HSD17B2','F-IER3','F-PTGDS','F-SFRP2'))
plot <- DotPlot(PURE, features = markers,cols = c("lightgrey","#FF0000"),col.min = 0,col.max = 2) + RotatedAxis() #c('white','#F8766D')
plot <- plot+scale_x_discrete(limits=markers) +labs(x=element_blank(), y=element_blank()) ### 调整cluster的顺序
plot[["guides"]][["size"]][["title"]] <- "Percent \nExpressed" ; plot[["guides"]][["colour"]][["title"]] <- "Average \nExpression"
ggsave("Mark_IFNs.png", plot = plot, width = 5.25, height = 4.75, dpi=300)


#######  作热图 趋化因子
markers <- c('CXCL14','CCL21','CCL19','CCL2','CXCL2','CXCL12',"CXCL16",'CXCL8','CCL11','CCL13')
col <- as.character(c('F-NGFR','F-IER3','F-ALPL','F-SFRP2','F-ISG15','F-STMN1','F-TPPP3','F-CPE','F-PCOLCE2',
                      'F-CXCL14','F-PTGDS','F-HSD17B2','F-POSTN')) 
Average <- AverageExpression(PURE, assays = NULL, features = NULL,
                             return.seurat = FALSE, add.ident = NULL, slot = "data",
                             use.scale = FALSE, use.counts = FALSE, verbose = TRUE)
gene0 <- Average$RNA %>% as.matrix()
TopMarkersExp <- gene0[rownames(gene0) %in% markers,]
TopMarkersExp <- TopMarkersExp[,col] # 调整热图的列，即调整数据框中cluster的排序
TopMarkersExp <- TopMarkersExp[markers,] # 调整热图的行，即调整数据框中markers的排序

plot <- pheatmap(TopMarkersExp,scale = "row", #clustering_method = "average",
                 cluster_rows = T,cluster_cols=F, # 对行和列进行聚类
                 border_color =NA, #"white", # 格子进行画线
                 cellwidth = 18, cellheight = 18,## 设置方格的大小
                 angle_col = 90,  ## 列的字体方向;没有对行的字体方向进行修改的
                 treeheight_row = 0, treeheight_col = 0, # 聚类树的高低，0为不画
                 fontsize = 10,fontsize_row = 12,fontsize_col = 14, # 字体大小
                 color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100)) # 调用颜色
ggsave("PURE_CKs.png", plot = plot, width = 10, height = 10, dpi=300)


### RNA Velocicty
source('/data2/teacher/yang/BCa/seurat/utils.R')
source('/data2/teacher/yang/BCa/seurat/seurat.R')
source('/data2/teacher/yang/BCa/seurat/velocyto.R')
library(ggplot2)
library(SeuratWrappers)

combined_looms <- readRDS("./combined_looms2.RDS")

dir.create('velocyto');setwd('velocyto')

set.seed(1200)
selected_so2 <- subset(PURE, downsample = 1200)

vel2 <- matrix_filter.Seurat(
  selected_so2, combined_looms,
  emat.min.max.cluster.average = 0.2, nmat.min.max.cluster.average = 0.05)

#### 计算RNA速率 velocyto2
rvel.cd2 <- gene.relative.velocity.estimates(
  # 输入对应细胞的emat、namt和cell.dist
  vel2$emat, vel2$nmat, cell.dist=vel2$cell.dist,
  deltaT=1, kCells=25, fit.quantile=0.02, n.cores=52)

#抽取出原始UMAP图中的对应颜色
colors <- DimPlot(object=selected_so2, reduction = "umap", label = TRUE) %>%
  extract_colors(vel2$emb)

# 绘图
show.velocity.on.embedding.cor(
  # 输入对应细胞的emb以及上一步计算出的RNA速率
  vel2$emb, rvel.cd2,
  n=2000, # 整体看箭头的趋势
  scale='sqrt', ###默认log,可选other values: 'sqrt','rank','linear'
  corr.sigma = 0.05, ##会影响箭头，默认0.05
  cex=0.8,##细胞圆点的大小
  arrow.scale=4,##箭头的长短
  show.grid.flow=T, ## FALSE是从每个细胞点向其他细胞指向
  cell.colors=ac(colors, alpha = 0.5),##  去颜色，可以显示黑白的
  min.grid.cell.mass=0.4, ##箭头数量，数值越大箭头越少
  fixed.arrow.length = F, ##T为固定显示箭头，默认为F
  plot.grid.points = F, ## T为背景网格点，默认为F
  grid.n=60,##  每个网格显示的箭头数量
  arrow.lwd=1, ##箭头加粗
  lwd=0,##每个细胞圆点阴影
  do.par=F,##参数不影响
  cell.border.alpha = 0.1, ##参数不影响
  n.cores=52 ) ##使用计数的CPU核心数 


###     SCENIC

library(SCENIC)
library(AUCell)
library(RcisTarget)
library(dplyr)
library(doRNG)
library(pheatmap)
#library(SCopeLoomR)
library(doMC)
library(rbokeh)
library(gistr)
library(limma)
library(tidyverse)
library(patchwork)

dir.create("SCENIC");setwd("SCENIC");dir.create("int")

SCENIC_test <- PURE

exprMat<-as.matrix(SCENIC_test@assays$RNA@counts)

cellInfo <- SCENIC_test@meta.data[,"new_seurat_clusters",drop=F]
head(cellInfo)
colnames(cellInfo)<- "CellType"
saveRDS(cellInfo, file="int/cellInfo.Rds")

n<-levels(cellInfo$CellType)%>%length()
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

for (i in 1:n){
  colVars <- list(CellType=c(i=gg_color_hue(i)))
}

names(colVars$CellType)<-as.character(levels(SCENIC_test@meta.data[["new_seurat_clusters"]]))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")

org <- "hgnc" 
dbDir <- './SCENIC/RcisTarget'  # RcisTarget databases location
myDatasetTitle <- "SCENIC_Test"     # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=48)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered <- log2(exprMat_filtered+1)
saveRDS(exprMat_filtered,file = 'int/exprMat_filtered.Rds')

runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions, nParts = 20)

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered,skipHeatmap=TRUE,skipTsne=TRUE)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_filtered)

##根据SCENIC结果作图
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled2 <- regulonActivity_byCellType_Scaled

####  删除extended的转录因子（该转录因子可靠性低）
regulonActivity_byCellType_Scaled2 <- regulonActivity_byCellType_Scaled2[-grep("_extended", rownames(regulonActivity_byCellType_Scaled2)),]

plot <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled2, #fontsize_row=12, 
                           color=colorRampPalette(c("navy","white","firebrick3"))(100), 
                           breaks=seq(-3, 3, length.out = 100),angle_col = 0,
                           treeheight_row=10, treeheight_col=10, border_color="white")
ggsave("output/SCENIC_1.png", plot = plot, width = 4.5, height = 15, dpi=300)  

dir.create('output/binaryAUC')

dr_coords <- Embeddings(SCENIC_test, reduction="umap") ## 基因的umap坐标
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[-grep("_extended",rownames(regulonActivity_byCellType_Scaled)),]
regulonNames <- rownames(regulonActivity_byCellType_Scaled)

dpi=300
for (i in 1:length(regulonNames)) {
  regulonNames[i]
  tf <- str_split(regulonNames[i], " ")[[1]][1]
  png(paste0("output/binaryAUC/binaryAUC_",tf,".png"), width=dpi*4.5,height=dpi*5,units = "px",res = dpi,type='cairo')
  AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tf), 
                  cex = 0.5, alphaOn = 1, alphaOff = 0.2,plots = "binaryAUC")
  dev.off()
  marker <- plot_density(PURE, reduction='umap', features = tf) 
  ggsave(paste0("output/binaryAUC/",tf,"_plot_density.png"), plot = marker, width = 4.5, height = 4, dpi=dpi)
}


##  GSVA
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(stringr)

dir.create('GSVA');setwd('GSVA')

scRNA <- PURE

# 选择基因集 HALLMARK
genesets <- msigdbr(species = "Homo sapiens", category = "H")
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol,genesets$gs_name)
expr <- AverageExpression(scRNA,assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]
expr <- as.matrix(expr)
gsva.res <- gsva(expr, genesets, method= "ssgsea")
rownames(gsva.res) <- gsub("HALLMARK_","",rownames(gsva.res))
plot <- pheatmap(gsva.res,show_colnames = T, scale = "row",
                 cluster_rows = TRUE, cluster_cols = TRUE,border_color = "white",angle_col = 90,
                 clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",
                 color = colorRampPalette(c("blue","white","red"))(100),#cutree_cols = 2,
                 #annotation_col = annotation_col,
                 cellwidth = 20, cellheight = 14, ## 设置方格的大小
                 fontsize_col = 14, # fontsize_row  行列的字体大小
                 treeheight_row = 5, treeheight_col = 5)
ggsave('Hallmark.png', plot = plot, width=10,height=20, dpi=300)

# 选择基因集 KEGG
genesets <- msigdbr(species = "Homo sapiens", category = "C2")
genesets <- subset(genesets, gs_subcat=="CP:KEGG",select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol,genesets$gs_name)
gsva.res2 <- gsva(expr, genesets, method= "ssgsea")
rownames(gsva.res2) <- gsub("KEGG_","",rownames(gsva.res2))
gsva.res2 <- gsva.res2[grepl("PATHWAY",rownames(gsva.res2)),]
rownames(gsva.res2) <- gsub("_SIGNALING_PATHWAY","",rownames(gsva.res2))
plot <- pheatmap(gsva.res2,show_colnames = T, scale = "row",
                 cluster_rows = TRUE, cluster_cols = TRUE,border_color = "white",angle_col = 90,
                 clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",
                 color = colorRampPalette(c("blue","white","red"))(100),#cutree_cols = 2,
                 #show_rownames=F,show_colnames=F,
                 cellwidth = 20, cellheight = 14, ## 设置方格的大小
                 fontsize_col = 14, # fontsize_row  行列的字体大小
                 treeheight_row = 5, treeheight_col = 5)
ggsave('KEGG2.png', plot = plot, width=10,height=10, dpi=300)

