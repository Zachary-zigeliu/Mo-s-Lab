.libPaths()
source('/data2/teacher/yang/R_package.R')

folder = "./Tissue"
setwd(folder);dir.create("RDS");dir.create("QC");dir.create("Script")
setwd("RDS")

object.list <- list(Tumor01,Core02,Core03,Edge02,Edge03,Tumor04,Tumor05,Tumor06,Tumor07,Tumor09,
                    Normal02,Normal03,Normal05,Normal06,Normal08,Normal09,Normal10,Normal11,Normal12)

# normalize and identify variable features for each dataset independently
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:50)
anchors
saveRDS(anchors,file = 'anchors.RDS')

Integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
Integrated
rm(object.list,anchors)

### first generate data and scale data in RNA assay
DefaultAssay(Integrated) <- "RNA"
Integrated[['percent.mt']] <- PercentageFeatureSet(Integrated, pattern = "^MT-")
Integrated <- NormalizeData(object = Integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
Integrated <- FindVariableFeatures(object = Integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
Integrated <- ScaleData(Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"))

##细胞周期相关基因比例
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
Integrated <- CellCycleScoring(Integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

## change to integrated assay
DefaultAssay(Integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
Integrated <- ScaleData(Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"))
Integrated <- RunPCA(Integrated, npcs = 100, verbose = FALSE)
plot <- ElbowPlot(Integrated, ndims = 100)


#样本信息
#合并样本-orig.ident_new2
Integrated[["orig.ident_new2"]]<-Integrated@meta.data$orig.ident
Integrated@meta.data$orig.ident_new2<-str_replace_all(Integrated@meta.data$orig.ident,
                                                      c("Normal10"="Patient10","Normal11"="Patient11", "Normal12"="Patient12", 
                                                        "Normal02"="Patient02", "Normal03"="Patient03", "Normal05"="Patient05", 
                                                        "Normal06"="Patient06", "Normal08"="Patient08", "Normal09"="Patient09",
                                                        "Tumor04"="Patient04", "Tumor05"="Patient05", 
                                                        "Tumor06"="Patient06", "Tumor07"="Patient07", "Tumor09"="Patient09",
                                                        "Tumor01"="Patient01", "Core02"="Patient02","Core03"="Patient03",
                                                        "Edge02"="Patient02","Edge03"="Patient03"))
Integrated@meta.data$orig.ident_new2 <- str_trim(Integrated@meta.data$orig.ident_new2,side = "right")
table(Integrated$orig.ident_new2)

Integrated[["orig.ident_new4"]]<-Integrated@meta.data$orig.ident
Integrated@meta.data$orig.ident_new4 <- str_replace_all(Integrated@meta.data$orig.ident,
                                                        c("Tumor01"="MIBC","Core02"="MIBC","Core03"="MIBC",
                                                          "Edge02"="MIBC","Edge03"="MIBC",
                                                          "Tumor04"="MIBC","Tumor05"="NMIBC","Tumor06"="NMIBC","Tumor07"="MIBC","Tumor09"="NMIBC",
                                                          "Normal02"="Normal","Normal03"="Normal","Normal05"="Normal",
                                                          "Normal06"="Normal","Normal08"="Normal","Normal09"="Normal",
                                                          "Normal10"="Normal", "Normal11"="Normal", "Normal12"="Normal"))
Integrated@meta.data$orig.ident_new4 <- str_trim(Integrated@meta.data$orig.ident_new4,side = "right")
table(Integrated@meta.data[["orig.ident_new2"]], Integrated$orig.ident_new4)
levels(Integrated)
table(Integrated$orig.ident_new4)

Idents(Integrated) <-Integrated@meta.data$orig.ident
saveRDS(Integrated,file = 'Tissue_Integrated.RDS')


setwd(folder);dir.create('QC');setwd('QC')

### umap
Integrated <- RunUMAP(Integrated, reduction = "pca", dims = 1:50)

#######  细胞大类
dir.create('CellTypes');setwd('CellTypes')

Tcells <-  c('CD2', 'CD3D','CD3E','CD3G') #darkorange 
Integrated$T.cell = apply(Integrated[['RNA']]@counts[c('CD2', 'CD3D','CD3E','CD3G'),],2, function(x) {ifelse(any(x>1),1,0)})
Mareker <- DimPlot(Integrated, reduction='umap', group.by='T.cell', cols=c('lightgrey','darkorange'))+ NoLegend()+ FigTheme2
ggsave("Tcells.png", plot = Mareker, width = 5, height = 5, dpi=300)

Bcells <-  c('CD79A','CD79B', 'CD19','MS4A1') #red 
Integrated$T.cell = apply(Integrated[['RNA']]@counts[c('CD79A','CD79B', 'CD19','MS4A1'),],2, function(x) {ifelse(any(x>1),1,0)})
Mareker <- DimPlot(Integrated, reduction='umap', group.by='T.cell', cols=c('lightgrey','red'))+ NoLegend()+ FigTheme2
ggsave("Bcells.png", plot = Mareker, width = 5, height = 5, dpi=300)

Myeloid <-  c('LYZ','CD68','CD14','S100A8','CPA3','CLC','MS4A6A','LAMP3','IL3RA')  #any(x>2) chocolate
Integrated$T.cell = apply(Integrated[['RNA']]@counts[c('LYZ','CD68','CD14','S100A8','CPA3','CLC','MS4A6A','LAMP3','IL3RA'),],2, function(x) {ifelse(any(x>2),1,0)})
Mareker <- DimPlot(Integrated, reduction='umap', group.by='T.cell', cols=c('lightgrey','chocolate'))+ NoLegend()+ FigTheme2
ggsave("Myeloid.png", plot = Mareker, width = 5, height = 5, dpi=300)

Fibro <-  c('COL1A1','COL4A1','COL3A1','DCN') #'FAP', any(x>2) blue
Integrated$T.cell = apply(Integrated[['RNA']]@counts[c('COL1A1','COL4A1','COL3A1','DCN'),],2, function(x) {ifelse(any(x>2),1,0)})
Mareker <- DimPlot(Integrated, reduction='umap', group.by='T.cell', cols=c('lightgrey','blue'))+ NoLegend()+ FigTheme2
ggsave("Fibro.png", plot = Mareker, width = 5, height = 5, dpi=300)

Pericyte <-  c('RGS5','MCAM','PDGFRB') #'ADRA2A', any(x>3) purple
Integrated$T.cell = apply(Integrated[['RNA']]@counts[c('RGS5','MCAM','PDGFRB'),],2, function(x) {ifelse(any(x>3),1,0)})
Mareker <- DimPlot(Integrated, reduction='umap', group.by='T.cell', cols=c('lightgrey','purple'))+ NoLegend()+ FigTheme2
ggsave("Pericyte.png", plot = Mareker, width = 5, height = 5, dpi=300)

Endothelium <-  c('FLT1','PLVAP','PROX1')  #deepskyblue
Integrated$T.cell = apply(Integrated[['RNA']]@counts[c('FLT1','PLVAP','PROX1'),],2, function(x) {ifelse(any(x>1),1,0)})
Mareker <- DimPlot(Integrated, reduction='umap', group.by='T.cell', cols=c('lightgrey','deepskyblue'))+ NoLegend()+ FigTheme2
ggsave("Endothelium.png", plot = Mareker, width = 5, height = 5, dpi=300)

SMC <-  c('ACTG2','DES','MYL9') #any(x>10) # limegreen
Integrated$T.cell = apply(Integrated[['RNA']]@counts[c('ACTG2','DES','MYL9'),],2, function(x) {ifelse(any(x>10),1,0)})
Mareker <- DimPlot(Integrated, reduction='umap', group.by='T.cell', cols=c('lightgrey','limegreen'))+ NoLegend()+ FigTheme2
ggsave("SMC.png", plot = Mareker, width = 5, height = 5, dpi=300)

Epithelium <-  c('EPCAM','KRT5', 'KRT17', 'KRT13','KRT20')  #darkgreen
Integrated$T.cell = apply(Integrated[['RNA']]@counts[c('EPCAM','KRT5', 'KRT17', 'KRT13','KRT20'),],2, function(x) {ifelse(any(x>1),1,0)})
Mareker <- DimPlot(Integrated, reduction='umap', group.by='T.cell', cols=c('lightgrey','darkgreen'))+ NoLegend()+ FigTheme2
ggsave("Epithelium.png", plot = Mareker, width = 5, height = 5, dpi=300)


#####  Samlpe origin
Idents(Integrated) <- Integrated@meta.data[["orig.ident_new4"]]
plot <- DimPlot(Integrated, reduction='umap',label = FALSE) +NoLegend() +FigTheme + FigAnno1 + FigAnno2
ggsave("Sample_origin.png", plot = plot, width = 9, height = 9, dpi=300)


### cluster
Integrated <- FindNeighbors(Integrated, reduction = "pca", dims = 1:50)
Integrated_1.4 <- FindClusters(Integrated, resolution = 1.4)

# find markers for every cluster compared to all remaining cells, report only the positive ones
Integrated_1.4.MAST <- FindAllMarkers(Integrated_1.4, assay = 'RNA', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST")


### 细胞比例图
RawCelltypes <- c('Fibroblasts',
                  'Epithelial cells','T cells','T cells','T cells','Fibroblasts',
                  'Epithelial cells','Mono-Macrophages','Epithelial cells','Mono-Macrophages','T cells',
                  'Epithelial cells','B cells','Neutrophils','Fibroblasts','T cells',
                  'Mono-Macrophages','Mast cells','T cells','T cells','Unassigned',
                  'T cells','Plasma cells','DCs','NK cells','Fibroblasts',
                  'Mono-Macrophages','SMCs','Endothelial cells','NK cells','Mono-Macrophages',
                  'T cells','T cells','Fibroblasts','Pericytes','T cells',
                  'T cells','Fibroblasts','DCs','DCs','DCs',
                  'Endothelial cells','Fibroblasts','T cells','Eosinophils','Neurons')
names(RawCelltypes) <- levels(Integrated_1.4)
Integrated_1.4 <- RenameIdents(Integrated_1.4, RawCelltypes)
id <- Integrated_1.4@active.ident
Integrated_1.4@meta.data$RawCelltypes <- id
Idents(Integrated_1.4) <- Integrated_1.4@meta.data[["seurat_clusters"]]

# CellAverage_1/CellProportion_1
Idents(Integrated_1.4) <- Integrated_1.4@meta.data[["RawCelltypes"]]
Idents(Integrated_1.4) <- factor(Idents(Integrated_1.4), 
                                 levels =c('T cells','NK cells','B cells','Plasma cells',
                                           'Mono-Macrophages','DCs','Mast cells','Neutrophils','Eosinophils',
                                           'Epithelial cells','Fibroblasts','Endothelial cells',
                                           'Pericytes','SMCs','Neurons','Unassigned'))

## Patient11的细胞数太少，后续将不纳入数据统计

######  细胞每样本平均数
Integrated_1.4 <- subset(Integrated_1.4, subset = orig.ident_new2 != 'Patient11') # Patient11细胞数太少，不纳入
Idents(Integrated_1.4) <- Integrated_1.4@meta.data[["RawCelltypes"]]

Normal <- subset(Integrated_1.4,orig.ident_new4=="Normal")
n <- levels(Normal@active.ident)
rep <- data.frame(n,rep("Normal",length(n)))
colnames(rep)<- c('Cluster','Type')
rep2 <- data.frame(table(Normal@active.ident)/8,stringsAsFactors = F) 
colnames(rep2) <- c('Cluster','Average')
Normal <- merge(rep2, rep)

MIBC <- subset(Integrated_1.4,orig.ident_new4=="MIBC")
n <- levels(MIBC@active.ident)
rep <- data.frame(n,rep("MIBC",length(n)))
colnames(rep)<- c('Cluster','Type')
rep2 <- data.frame(table(MIBC@active.ident)/7,stringsAsFactors = F) 
colnames(rep2) <- c('Cluster','Average')
MIBC <- merge(rep2, rep)

NMIBC <- subset(Integrated_1.4,orig.ident_new4=="NMIBC")
n <- levels(NMIBC@active.ident)
rep <- data.frame(n,rep("NMIBC",length(n)))
colnames(rep)<- c('Cluster','Type')
rep2 <- data.frame(table(NMIBC@active.ident)/3,stringsAsFactors = F) 
colnames(rep2) <- c('Cluster','Average')
NMIBC <- merge(rep2, rep)

Ave2 <- rbind(Normal,MIBC) %>% rbind(NMIBC)  

library(RColorBrewer)
nb.cols <- length(levels(Integrated_1.4@active.ident))
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
  scale_x_discrete(limits=c("Normal", "MIBC","NMIBC"))+
  scale_fill_manual(values = mycolors)+ ##调用颜色方案
  RotatedAxis()
ggsave("CellAverage_B.png", plot = plot, width = 5, height = 6, dpi=300)  


######  细胞比例
cell.prop<-as.data.frame(prop.table(table(Idents(Integrated_1.4), Integrated_1.4$orig.ident_new4)))
colnames(cell.prop)<-c("Cluster","origin","Proportion")

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
  scale_x_discrete(limits=c("Normal", "MIBC","NMIBC"))+
  scale_y_continuous(labels = scales::percent)+ #纵坐标变为百分比
  scale_fill_manual(values = mycolors) + ##调用颜色方案
  RotatedAxis()
ggsave("CellProportion_B.png", plot = plot, width = 5, height = 6, dpi=300)  

## Patient Proportion  
Normal <- subset(Integrated_1.4, orig.ident_new4 == "Normal")
cell.prop_Normal<-as.data.frame(prop.table(table(Idents(Normal), Normal$orig.ident_new2)))
colnames(cell.prop_Normal)<-c("Celltype","origin","Proportion")
cell.prop_Normal$Tissue <- "Normal"
plot <- ggplot(cell.prop_Normal,aes(origin,Proportion,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  theme(panel.background = element_rect(fill = NA,colour = 'black'),axis.ticks = element_blank(),#axis.line = element_line(color = 'black'),
        legend.position='none',
        axis.title = element_blank(),
        axis.text.x =element_text(size = 14, color = 'black'),axis.text.y =element_text(size = 12, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 20))+
  scale_y_continuous(labels = scales::percent)+ #纵坐标变为百分比
  scale_fill_manual(values = mycolors) + ##调用颜色方案
  RotatedAxis()
ggsave("PatientProp_Normal.png", plot = plot, width = 5.25, height = 5.5, dpi=300)  

MIBC <- subset(Integrated_1.4, orig.ident_new4 == "MIBC")
cell.prop_MIBC<-as.data.frame(prop.table(table(Idents(MIBC), MIBC$orig.ident_new2)))
colnames(cell.prop_MIBC)<-c("Celltype","origin","Proportion")
cell.prop_MIBC$Tissue <- "MIBC"
plot <- ggplot(cell.prop_MIBC,aes(origin,Proportion,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  theme(panel.background = element_rect(fill = NA,colour = 'black'),axis.ticks = element_blank(),#axis.line = element_line(color = 'black'),
        legend.position='none',
        axis.title = element_blank(),
        axis.text.x =element_text(size = 14, color = 'black'),axis.text.y =element_text(size = 12, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 20))+
  scale_y_continuous(labels = scales::percent)+ #纵坐标变为百分比
  scale_fill_manual(values = mycolors) + ##调用颜色方案
  RotatedAxis()
ggsave("PatientProp_MIBC.png", plot = plot, width = 3.6, height = 5.5, dpi=300)  

NMIBC <- subset(Integrated_1.4, orig.ident_new4 == "NMIBC")
cell.prop_NMIBC<-as.data.frame(prop.table(table(Idents(NMIBC), NMIBC$orig.ident_new2)))
colnames(cell.prop_NMIBC)<-c("Celltype","origin","Proportion")
cell.prop_NMIBC$Tissue <- "NMIBC"
plot <- ggplot(cell.prop_NMIBC,aes(origin,Proportion,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  theme(panel.background = element_rect(fill = NA,colour = 'black'),axis.ticks = element_blank(),#axis.line = element_line(color = 'black'),
        legend.position='none',
        axis.title = element_blank(),
        axis.text.x =element_text(size = 14, color = 'black'),axis.text.y =element_text(size = 12, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 20))+
  scale_y_continuous(labels = scales::percent)+ #纵坐标变为百分比
  scale_fill_manual(values = mycolors) + ##调用颜色方案
  RotatedAxis()
ggsave("PatientProp_NMIBC.png", plot = plot, width = 2.5, height = 5.5, dpi=300)  


###  -------------细胞比例箱线图-----------
Idents(Integrated_1.4) <- Integrated_1.4@meta.data$seurat_clusters
RawCelltypes2 <- c('Fibroblasts',
                   'Epithelial cells','Immune cells','Immune cells','Immune cells','Fibroblasts',
                   'Epithelial cells','Immune cells','Epithelial cells','Immune cells','Immune cells',
                   'Epithelial cells','Immune cells','Immune cells','Fibroblasts','Immune cells',
                   'Immune cells','Immune cells','Immune cells','Immune cells','Unassigned',
                   'Immune cells','Immune cells','Immune cells','Immune cells','Fibroblasts',
                   'Immune cells','SMCs','Endothelial cells','Immune cells','Immune cells',
                   'Immune cells','Immune cells','Fibroblasts','Pericytes','Immune cells',
                   'Immune cells','Fibroblasts','Immune cells','Immune cells','Immune cells',
                   'Endothelial cells','Fibroblasts','Immune cells','Immune cells','Neurons')
names(RawCelltypes2) <- levels(Integrated_1.4)
Integrated_1.4 <- RenameIdents(Integrated_1.4, RawCelltypes2)
id <- Integrated_1.4@active.ident
Integrated_1.4@meta.data$RawCelltypes2 <- id
Idents(Integrated_1.4) <- Integrated_1.4@meta.data[["seurat_clusters"]]

Integrated_1.4 <- subset(Integrated_1.4, subset = orig.ident_new2 != 'Patient11')
Idents(Integrated_1.4) <- Integrated_1.4@meta.data[["RawCelltypes2"]]
Idents(Integrated_1.4) <- factor(Idents(Integrated_1.4), 
                                 levels =c("Immune cells",
                                           "Epithelial cells","Fibroblasts","Endothelial cells","Pericytes","SMCs",
                                           "Neurons",'Unassigned'))

Cell.Num <- table(Idents(Integrated_1.4), Integrated_1.4$orig.ident_new2, Integrated_1.4$orig.ident_new4)
data.m <- melt(Cell.Num)
colnames(data.m)<-c("Celltype","Patient","origin","CellNum")

Patients <- c('Patient01','Patient02','Patient03','Patient04','Patient05',
              'Patient06','Patient07','Patient08','Patient09','Patient10','Patient12')
origins <- c('NMIBC','MIBC','Normal')
P_T <- table(Integrated_1.4$orig.ident_new2,Integrated_1.4$orig.ident_new4)

### 占各样品细胞总数的比例
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
Subtypes <- as.character(levels(Integrated_1.4))

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
  plot <- g + geom_signif(comparisons = compared_list, step_increase = 0.2) 
  ggsave(paste0(Subtype,'_boxplot_A.png'), plot = plot, width = 3.25, height = 4.5, dpi=300)  
} 

ImmuneCell <- data.M[which(data.M$Celltype %in% "Immune cells"),]

g <- ggplot(data=ImmuneCell,aes(x=origin,y=Prop,colour=origin)) + 
  stat_boxplot(geom = "errorbar",width=0.15)+
  labs(x = element_blank(), y = 'Cell proportion', title = 'Immune cells')+
  geom_boxplot() + 
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),  # panel.background = element_rect(fill = NA,colour = 'black'),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.title =element_text(hjust = 0.5,size = 16),
        axis.text.x =element_text(size = 14, color = 'black'),axis.text.y =element_text(size = 14, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 18)) + NoLegend()+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(limits=c("Normal", "MIBC","NMIBC"))
plot <- g + geom_signif(comparisons = compared_list, step_increase = 0.2)
ggsave("ImmueCell_boxplot_Tissue_A.png", plot = plot, width = 3.5, height = 4, dpi=300)  







































