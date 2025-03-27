source('/data2/teacher/yang/R_package.R')

########  Tcell  ########
Cell <- "Tcell"
dir = "~/Tissue/re1.4"  ## 新服务器路径
folder <- paste0(dir,"/",Cell)
res <- "re0.8"

Integrated_1.4 <- readRDS("~/Tissue/RDS/Tissue_Integrated_1.4.RDS")

Sub <- subset(Integrated_1.4, idents = c(2,3,4,10,15,18,19,21,31,32,35,36,43))  ## 提取总图中相应的细胞类群

##重新分类
Sub <- NormalizeData(Sub, normalization.method = "LogNormalize", scale.factor = 10000)
Sub <- FindVariableFeatures(Sub, selection.method = "vst", nfeatures = 2000)
Sub <- ScaleData(Sub, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE ) # "percent.mt","percent.ribo"
Sub <- RunPCA(Sub, npcs = 100, verbose = FALSE)

## Harmony去批次效应
Sub <- Sub %>% RunHarmony("orig.ident", plot_convergence = TRUE)

Sub <- RunUMAP(Sub, reduction = "harmony", dims = 1:30)
Sub <- FindNeighbors(Sub, reduction = "harmony", dims = 1:30)
Sub_0.8 <- FindClusters(Sub, resolution = 0.8)
Sub_0.8.MAST <- FindAllMarkers(Sub_0.8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST")

## 去除混杂细胞后，提取相应的细胞类群
clus <- c(0,1,2,3,4,5,6,8,9,10)
PURE <- subset(Sub_0.8, idents = clus) 
PURE@meta.data[["seurat_clusters"]] <- droplevels(PURE@meta.data[["seurat_clusters"]])

## 命名 new_seurat_clusters
new.cluster.ids <- c('0','1','2','3','4','5','6','7','8','9')
names(new.cluster.ids) <- levels(PURE)
PURE <- RenameIdents(PURE, new.cluster.ids)
PURE@meta.data$new_seurat_clusters <- PURE@active.ident
PURE <- RunUMAP(PURE, reduction = "harmony", dims = 1:30)

## 命名 new_seurat_clusters2
new.cluster.ids <- c('CD8-GZMK','CD8-GNLY','CD4-FOXP3','Naive-IL7R',
                     'CD8-HAVCR2','CD8-HSPA1A','CD4-CXCL13','CD8-MKI67',
                     'CD8-STMN1','CD8-ISG15')
names(new.cluster.ids) <- levels(PURE)
PURE <- RenameIdents(PURE, new.cluster.ids)
PURE@meta.data$new_seurat_clusters2 <- PURE@active.ident

## 命名 new_seurat_clusters3
new.cluster.ids <- c('0    CD8-GZMK','1    CD8-GNLY','2    CD4-FOXP3','3    Naive-IL7R',
                     '4    CD8-HAVCR2','5    CD8-HSPA1A','6    CD4-CXCL13','7    CD8-MKI67',
                     '8    CD8-STMN1','9    CD8-ISG15')
names(new.cluster.ids) <- levels(PURE)
PURE <- RenameIdents(PURE, new.cluster.ids)
PURE@meta.data$new_seurat_clusters3 <- PURE@active.ident

plot <- DimPlot(PURE, reduction = "umap", label = TRUE, label.size = 6)
ggsave("PURE_new3A.png", plot = plot, width = 7.5, height = 6, dpi=300)

Idents(PURE) <- PURE@meta.data[["new_seurat_clusters"]]
saveRDS(PURE,file = paste0(folder,"/RDS/PURE.RDS"))


## 热图 细胞cluster特异marker
Markers <- c("CD4","CD8A","CD8B",
             "SELL","LEF1","TCF7","IL7R","CCR7", # Naive
             "HAVCR2","LAG3","TIGIT","CTLA4","PDCD1", # Tex 共抑制分子
             "GZMK","GZMA","NKG7","GZMB","GNLY","PRF1", # CTL 
             "IL2RA","FOXP3","IKZF2", # Treg
             "PCNA","STMN1","MKI67","UBE2C") # Cycling

Idents(PURE) <- PURE@meta.data[["new_seurat_clusters2"]]
Average <- AverageExpression(PURE, assays = NULL, features = NULL,
                             return.seurat = FALSE, add.ident = NULL, slot = "data",
                             use.scale = FALSE, use.counts = FALSE, verbose = TRUE)
gene1 <- Average$RNA %>% as.matrix()
TopMarkersExp <- gene1[rownames(gene1) %in% Markers, ]

col <- c('Naive-IL7R','CD4-FOXP3','CD4-CXCL13','CD8-ISG15','CD8-GNLY','CD8-HSPA1A','CD8-HAVCR2',
         'CD8-STMN1','CD8-MKI67','CD8-GZMK') 
TopMarkersExp <- TopMarkersExp[,col] # 调整热图的列，即调整数据框中cluster的排序
TopMarkersExp <- TopMarkersExp[Markers,] # 调整热图的行，即调整数据框中Markers的排序
plot <- pheatmap(TopMarkersExp,scale = "row", 
                 cluster_rows = F,cluster_cols=F, # 对行和列进行聚类
                 border_color =NA, #"white", # 格子进行画线
                 gaps_row = c(3,8,13,19,22),gaps_col = c(1,3,7,9),
                 angle_col = 90,  ## 列的字体方向;没有对行的字体方向进行修改的
                 fontsize = 10,fontsize_row = 10,fontsize_col = 10, # 字体大小
                 color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100)) # 调用颜色
ggsave("PURE_markgenes2.png", plot = plot, width = 4.5, height = 5.5, dpi=300)


######  细胞每样本平均数
PURE <- subset(PURE, subset = orig.ident_new2 != 'Patient11') # 不纳入Patient11，少一个Normal样本

Normal <- subset(PURE,orig.ident_new4=="Normal")
n <- levels(Normal@active.ident)
rep <- data.frame(n,rep("Normal",length(n)))
colnames(rep)<- c('Cluster','Type')
rep2 <- data.frame(table(Normal@active.ident)/8,stringsAsFactors = F) # 不纳入Patient11，少一个Normal样本 
colnames(rep2) <- c('Cluster','Average')
Normal <- merge(rep2, rep)

MIBC <- subset(PURE,orig.ident_new4=="MIBC")
n <- levels(MIBC@active.ident)
rep <- data.frame(n,rep("MIBC",length(n)))
colnames(rep)<- c('Cluster','Type')
rep2 <- data.frame(table(MIBC@active.ident)/7,stringsAsFactors = F) 
colnames(rep2) <- c('Cluster','Average')
MIBC <- merge(rep2, rep)

NMIBC <- subset(PURE,orig.ident_new4=="NMIBC")
n <- levels(NMIBC@active.ident)
rep <- data.frame(n,rep("NMIBC",length(n)))
colnames(rep)<- c('Cluster','Type')
rep2 <- data.frame(table(NMIBC@active.ident)/3,stringsAsFactors = F) 
colnames(rep2) <- c('Cluster','Average')
NMIBC <- merge(rep2, rep)

Ave2 <- rbind(Normal,MIBC) %>% rbind(NMIBC) 

nb.cols <- length(levels(PURE@active.ident))
mycolors <- colorRampPalette(brewer.pal(nb.cols, "Set3"))(nb.cols)

#theme_bw()加框
plot <- ggplot(Ave2,aes(Type,Average,fill=Cluster))+
  geom_bar(stat='identity', position='stack')+# 'dodge'为并排式
  #theme_bw() +
  labs(x = '',y = 'Average number (Cells/Sample)', title = '')+
  theme(panel.background = element_rect(fill = NA,colour = 'black'),#"transparent"
        axis.ticks = element_blank(),#axis.line = element_blank(),
        legend.text = element_text(size=14),
        axis.title =element_text(size = 20),
        axis.text.x =element_text(size = 14, color = 'black'),axis.text.y =element_text(size = 12, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 20))+
  scale_fill_manual(values = mycolors)+ ##调用颜色方案
  RotatedAxis()#+
# guides(fill=guide_legend(title=NULL))
ggsave("CellAverage2.png", plot = plot, width = 4.5, height = 6, dpi=300)  


### 各亚群的细胞比例
Cell.Num <- table(Idents(PURE), PURE$orig.ident_new2, PURE$orig.ident_new4)
data.m <- melt(Cell.Num)
colnames(data.m)<-c("Celltype","Patient","origin","CellNum")

Patients <- c('Patient01','Patient02','Patient03','Patient04','Patient05',
              'Patient06','Patient07','Patient08','Patient09','Patient10','Patient12')
origins <- c('NMIBC','MIBC','Normal')
P_T <- table(Integrated_1.4$orig.ident_new2,Integrated_1.4$orig.ident_new4)

nb.cols <- length(levels(PURE@active.ident))
mycolors <- colorRampPalette(brewer.pal(nb.cols, "Set3"))(nb.cols)

data.M <- NULL
for (p in Patients) {
  data.m2 <- data.m[which(data.m$Patient %in% p),]
  for(origin in origins){
    data.m3 <-  data.m2[which(data.m2$origin %in% origin),]
    data.m3$Prop <- data.m3$CellNum
    data.m3$Prop[1:nb.cols] <- data.m3$Prop[1:nb.cols]/as.numeric(P_T[p,origin]) ## 注意细胞类群的数据
    data.M <- rbind(data.M,data.m3)
  }
}

data.M <- data.M[-which(data.M$Prop %in% NaN),]
data.M$CellNum <- NULL

compared_list <- list(c("Normal", "MIBC"), c("MIBC", "NMIBC"), c("Normal", "NMIBC"))
Subtypes <- as.character(levels(PURE))

library(ggsignif)#加载显著性包

for (i in 1:nb.cols) {
  Subtype <- Subtypes[i]
  data_plot <- data.M[which(data.M$Celltype %in% Subtype),]
  g <- ggplot(data_plot,aes(x=origin,y=Prop,colour=origin)) + 
    stat_boxplot(geom = "errorbar",width=0.15)+
    labs(x = element_blank(), y = 'Cell proportion',title = Subtype)+ geom_boxplot() + 
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
          panel.background = element_blank(),  # panel.background = element_rect(fill = NA,colour = 'black'),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_blank(), axis.title =element_text(hjust = 0.5,size = 16),
          axis.text =element_text(size = 14, color = 'black'),
          plot.title =element_text(hjust = 0.5, size = 16)) + NoLegend()+
    scale_y_continuous(labels = scales::percent)+RotatedAxis()+
    scale_x_discrete(limits=c("Normal", "MIBC","NMIBC"))
  plot <- g + geom_signif(comparisons = compared_list, step_increase = 0.17) #"anova" test = kruskal.test,
  ggsave(paste0(Subtype,'_boxplot_A.png'), plot = plot, width = 3.25, height = 4.5, dpi=300)  
} 


## 气泡图 （细胞因子）
markers <- c("CCL3","CCL4","CCL4L2","CCL5","CXCL13","IFNG","CSF1","TGFB1",
             "CCR4","CCR5","CCR6","CCR7","CXCR3","CXCR4","CXCR6","IFNGR1","TGFBR2",
             "ICAM1","ICAM2","ICAM3","SELL",
             "CD69","ITGA1","ITGA2","ITGA4","ITGAE","ITGAL","ITGB1","ITGB2","ITGB7")
plot <- DotPlot(PURE, features = markers,cols = c("lightgrey","#FF0000"),col.min = 0,col.max = 2) + RotatedAxis() #c('white','#F8766D')
plot <- plot+scale_x_discrete(limits=markers) +labs(x=element_blank(), y=element_blank()) ### 调整cluster的顺序
plot
ggsave("Mark_CK.png", plot = plot, width = 12, height = 4, dpi=300)


##  小提琴图   
Idents(PURE) <- PURE@meta.data[["orig.ident_new4"]]
Idents(PURE) <- factor(Idents(PURE), levels =c('Normal','MIBC','NMIBC'))
markers <- c("PDCD1","BTLA","CD244","CEACAM1","HAVCR2","LAG3","CTLA4","TIGIT", ## 抑制型共刺激分子
             "TNFRSF9","TNFRSF4", "TNFRSF18","ICOS","CD27","CD28")  ## 激活型共刺激分子
for(marker in markers){
  plot <- VlnPlot(PURE, features = marker,pt.size = 0)+NoLegend()+
    xlab(element_blank()) + ylab(marker) + ggtitle(element_blank()) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),#legend.position = "none",
          axis.text.y = element_blank(),axis.title.y =element_blank())
  ggsave(paste0("violin_",marker,".png"), plot = plot, width = 2, height = 0.85, dpi=300)
}

##  小提琴图   在MIBC和NMIBC中比较各细胞类群基因的表达 
Idents(PURE) <- PURE@meta.data[["new_seurat_clusters2"]]
Tumor <- subset(PURE, subset = orig.ident_new == "Tumor")
##调整数据中细胞类群顺序
Idents(Tumor) <- factor(Idents(Tumor), 
                        levels =c('CD4-FOXP3','CD4-CXCL13','CD8-MKI67','CD8-STMN1','CD8-HAVCR2','CD8-HSPA1A','CD8-ISG15',
                                  'CD8-GNLY','CD8-GZMK','Naive-IL7R'))
markers <- c("PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","BTLA","CD244","CEACAM1", ## 抑制型共刺激分子
             "TNFRSF9","TNFRSF18","TNFRSF4", "ICOS","CD27","CD28")  ## 激活型共刺激分子
for(marker in markers){
  plot <- VlnPlot(Tumor, features = marker,split.by = "orig.ident_new4",pt.size = 0)+NoLegend()+
    xlab(element_blank()) + ylab(marker) + ggtitle(element_blank()) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(),axis.title.y =element_blank())
  ggsave(paste0("violin_",marker,".png"), plot = plot, width = 4.5, height = 1, dpi=300)
}


##### RNA Velocicty
source('/data2/teacher/yang/BCa/seurat/utils.R')
source('/data2/teacher/yang/BCa/seurat/seurat.R')
source('/data2/teacher/yang/BCa/seurat/velocyto.R')
library(ggplot2)
library(SeuratWrappers)

combined_looms <- readRDS("./combined_looms2.RDS")

dir.create('velocyto');setwd('velocyto')

set.seed(1500)
selected_so2 <- subset(PURE, downsample = 1500) 

vel2 <- matrix_filter.Seurat(
  selected_so2, combined_looms,
  emat.min.max.cluster.average = 0.2, nmat.min.max.cluster.average = 0.05)

rvel.cd2 <- gene.relative.velocity.estimates(
  vel2$emat, vel2$nmat, cell.dist=vel2$cell.dist,
  deltaT=1, kCells=25, fit.quantile=0.02, n.cores=1)

colors <- DimPlot(object=selected_so2, reduction = "umap", label = TRUE) %>%
  extract_colors(vel2$emb)

# 绘图
p2 <- show.velocity.on.embedding.cor(
  vel2$emb, rvel.cd2,
  n=2000, # 整体看箭头的趋势
  scale='sqrt', ###默认log,可选other values: 'sqrt','rank','linear'
  corr.sigma = 0.05, ##会影响箭头，默认0.05
  cex=0.8,##细胞圆点的大小
  arrow.scale=6,##箭头的长短
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
  n.cores=48 ) ##使用计数的CPU核心数 


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






