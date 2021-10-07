#####################################
# Author: Wilson Tan
# Last Updated: 7 Oct 2021
#####################################  

## Dependencies / packages
require(Seurat)
require(patchwork)

## Reading in datasets from individual analysis for integration
sham_p1 = readRDS("../JXZ_30_CKDL210017820-1a-SI_TT_A5_HH5T2DSX2/majorHEART.rds")
sham_p1@meta.data$condition = rep("SHAM_P1")

mi_p1 = readRDS("../JXZ_31_CKDL210017821-1a-SI_TT_A6_HH5T2DSX2//majorHEART.rds")
mi_p1@meta.data$condition = rep("MI_P1")

sham_p9 = readRDS("../JXZ_37_CKDL210017826-1a-SI_TT_A7_HH5T2DSX2/majorHEART.rds")
sham_p9@meta.data$condition = rep("SHAM_P9")

mi_p9 = readRDS("../JXZ_38_CKDL210017827-1a-SI_TT_A8_HH5T2DSX2/majorHEART.rds")
mi_p9@meta.data$condition = rep("MI_P9")


## Form a list with object that you intend to merge
heart.list <- list(sham_p1=sham_p1,mi_p1=mi_p1,sham_p9=sham_p9,mi_p9=mi_p9)
heart.list

## Identify top 2000 most variable features
features <- SelectIntegrationFeatures(object.list=heart.list)
features

## Form an anchor
heart.anchors <- FindIntegrationAnchors(object.list=heart.list, anchor.features=features)
save(heart.anchors,features,file="heart_step1.RData")


## Data integration
load(file="heart_step1.RData")
heart.combined <- IntegrateData(anchorset=heart.anchors)
save(heart.combined,features,file="heart_step2.RData")

## Rerun classification
load("heart_step2.RData")
DefaultAssay(heart.combined) <- "integrated"
heart.combined <- ScaleData(heart.combined)
heart.combined <- RunPCA(heart.combined,npcs=100)
heart.combined <- RunUMAP(heart.combined,reduction="pca",dims=1:100)
heart.combined <- FindNeighbors(heart.combined,reduction="pca",dims=1:100)
heart.combined <- FindClusters(heart.combined,resolution=1)
save(heart.combined,file="heart_step3.RData")


## run UMAP plotting
load("heart_step3.RData")
options(bitmapType="cairo")
png("UMAP_condition.png",width=3000,height=1500,res=300)
p1 <- DimPlot(heart.combined, reduction = "umap", group.by = "condition")
p2 <- DimPlot(heart.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
dev.off()

png("UMAP_split_condition.png",width=3000,height=3000,res=300)
DimPlot(heart.combined, reduction = "umap", split.by = "condition",ncol=2)
dev.off()


## Plot features based on defined markers
## To do: create a list of markers to cell type table
load("heart_step3.RData")
png("MARKER_featureplot.png",width=6000,height=6000,res=300)
FeaturePlot(heart.combined, features=c("Myl2","Tnnt2","Myh6","Nppa","Dach1","Emcn","Egfl7","Vwf","Cdh5","Tie1","Postn","Dcn","Fn1","Rgs5","Abcc9","Pcdh7","Cd163","Mrc1","Ikzf1","Fyb1"),raster=FALSE)
dev.off()

## Once you are happy with the plotting
## manually assign the cells to known type
## if unknown, leave it as OTHERS
heart <- RenameIdents(
        object = heart.combined,
        "0"="CM",
        "1"="CM",
        "2"="CM",
        "3"="FB",
        "4"="ENDO",
        "5"="CM",
        "6"="MACRO",
        "7"="CM",
        "8"="dCM",
        "9"="FB",
        "10"="PERI",
        "11"="FB",
        "12"="OTHER1",
        "13"="FB",
        "14"="FB",
        "15"="CM",
        "16"="ENDO",
        "17"="ENDO",
        "18"="ENDO",
	"19"="CM",
	"20"="PERI",
	"21"="OTHER2",
	"22"="PERI",
	"23"="ENDO",
	"24"="OTHER3",
	"25"="ENDO"
)

### Basically replot
umap.coord = as.data.frame(heart@reductions$umap@cell.embeddings)
umap.coord$group = as.factor(Idents(heart))
umap.coord$condition = heart$condition
umap.coord$condition = factor(umap.coord$condition,levels=c("SHAM_P1","MI_P1","SHAM_P9","MI_P9"))

options(bitmapType="cairo")
png("UMAP_namedTYPE.png",width=3000,height=3000,res=300)
p1 <- ggplot(umap.coord,aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color=group),size=1.2)+theme_bw()+theme(panel.grid=element_blank())
p1 <- p1 + scale_color_manual(values=c("darkred","red","orange","yellow","turquoise","blue","pink","black","brown","purple","pink","salmon","bisque","grey30","grey60","green","darkgreen","yellowgreen"))
print(p1)
dev.off()
png("UMAP_namedTYPE_byCondition.png",width=6000,height=3000,res=300)
p1 <- ggplot(umap.coord,aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color=group),size=1.2)+theme_bw()+theme(panel.grid=element_blank())
p1 <- p1 + facet_grid(.~condition,scale="free_y")
p1 <- p1 + scale_color_manual(values=c("darkred","red","orange","yellow","turquoise","blue","pink","black","brown","purple","pink","salmon","bisque","grey30","grey60","green","darkgreen","yellowgreen"))
print(p1)
dev.off()


saveRDS(heart,file="defined.celltype.RDS")


### THIS PATRT IS DIFFERENTIAL EX ANALYSIS
## Comparing MI vs Sham for P1
DefaultAssay(heart) <- "RNA"
heart$celltype.MI <- paste(Idents(heart), heart$condition, sep = "_")
heart$celltype <- Idents(heart)
Idents(heart) <- "celltype.MI"

group.list = unique(gsub("_MI_P[0-9]","",gsub("_SHAM_P[0-9]","",Idents(heart))))
group.list

for(i in 1:length(group.list))
{
        print(group.list[i])
        mi.response <- FindMarkers(heart, ident.1 = paste(group.list[i],"_SHAM_P1",sep=""), ident.2 = paste(group.list[i],"_MI_P1",sep=""), verbose = FALSE)
        write.table(mi.response,file=paste(group.list[i],"__P1_DE.txt",sep=""),sep="\t",quote=FALSE)
        print(table(mi.response$p_val_adj<0.05))


}

for(i in 1:length(group.list))
{
        print(group.list[i])
        mi.response <- FindMarkers(heart, ident.1 = paste(group.list[i],"_SHAM_P9",sep=""), ident.2 = paste(group.list[i],"_MI_P9",sep=""), verbose = FALSE)
        write.table(mi.response,file=paste(group.list[i],"__P9_DE.txt",sep=""),sep="\t",quote=FALSE)
        print(table(mi.response$p_val_adj<0.05))


}

#### THIS PART IS TRYING TO SEE IF CELL PROPORTION CHANGE
### Count proportion of cell type
require(stringr)
require(ggplot2)
options(bitmapType="cairo")
celltype.stats = as.data.frame(table(Idents(heart)))
celltype.stats = cbind(celltype.stats,as.data.frame(matrix(unlist(str_split(celltype.stats$Var1,"_")),byrow=TRUE,ncol=3)))
names(celltype.stats)=c("ID","Freq","CT","Condition","Age")
celltype.stats$id = paste(celltype.stats$Condition,celltype.stats$Age)

total = aggregate(celltype.stats$Freq,by=list(celltype.stats$Condition,celltype.stats$Age),FUN=sum)
names(total) = c("Condition","Age","Total")
total$id = paste(total$Condition,total$Age)

celltype.stats$Total = total$Total[match(celltype.stats$id,total$id)]
celltype.stats$percent = celltype.stats$Freq/celltype.stats$Total*100

celltype.stats$Condition = factor(celltype.stats$Condition,levels=c("SHAM","MI"))
png("celltype.proportion.byCT.png",width=1500,height=1500,res=300)
p1 <- ggplot(celltype.stats,aes(x=Condition,y=percent))+geom_bar(stat="identity",position="stack",aes(fill=CT),width=0.75)+theme_bw()+facet_grid(.~Age,space="free",scale="free")
p1 <- p1 + theme(panel.grid=element_blank(),legend.position="bottom") + scale_fill_manual(values=c("darkred","red","orange","yellow","turquoise","blue","pink","black","purple","purple","pink","salmon","bisque","grey30","grey60","green","darkgreen","yellowgreen"))
p1 <- p1 + xlab("") + ylab("Percentage")
print(p1)
dev.off()


write.table(celltype.stats,file="celltype.stats.txt",sep="\t",quote=FALSE,row.names=F)



#### THIS PART RUN DE ANALYSIS ON EVERY SUBTYPE
## Comparing MI vs Sham for P1 per subtype
heart = readRDS("defined.celltype.RDS")
DefaultAssay(heart) <- "RNA"
heart$celltype.MI <- paste(Idents(heart), heart$integrated_snn_res.1,heart$condition,sep="_")
heart$celltype <- Idents(heart)
Idents(heart) <- "celltype.MI"

group.list = unique(gsub("_MI_P[0-9]","",gsub("_SHAM_P[0-9]","",Idents(heart))))
group.list

for(i in 1:length(group.list))
{
        print(group.list[i])
        mi.response <- FindMarkers(heart, ident.1 = paste(group.list[i],"_SHAM_P1",sep=""), ident.2 = paste(group.list[i],"_MI_P1",sep=""), verbose = FALSE)
        write.table(mi.response,file=paste(group.list[i],"__P1_DE.txt",sep=""),sep="\t",quote=FALSE)
        print(table(mi.response$p_val_adj<0.05))


}

for(i in 1:length(group.list))
{
        print(group.list[i])
        mi.response <- FindMarkers(heart, ident.1 = paste(group.list[i],"_SHAM_P9",sep=""), ident.2 = paste(group.list[i],"_MI_P9",sep=""), verbose = FALSE)
        write.table(mi.response,file=paste(group.list[i],"__P9_DE.txt",sep=""),sep="\t",quote=FALSE)
        print(table(mi.response$p_val_adj<0.05))


}

### Count proportion of cell type
require(stringr)
require(ggplot2)
options(bitmapType="cairo")
celltype.stats = as.data.frame(table(Idents(heart)))
celltype.stats = cbind(celltype.stats,as.data.frame(matrix(unlist(str_split(celltype.stats$Var1,"_")),byrow=TRUE,ncol=4)))
names(celltype.stats)=c("ID","Freq","CT","ST","Condition","Age")
celltype.stats$id = paste(celltype.stats$Condition,celltype.stats$Age)

total = aggregate(celltype.stats$Freq,by=list(celltype.stats$Condition,celltype.stats$Age),FUN=sum)
names(total) = c("Condition","Age","Total")
total$id = paste(total$Condition,total$Age)

celltype.stats$Total = total$Total[match(celltype.stats$id,total$id)]
celltype.stats$percent = celltype.stats$Freq/celltype.stats$Total*100

celltype.stats$CTST = paste(celltype.stats$CT,celltype.stats$ST,sep="_")

celltype.stats$Condition = factor(celltype.stats$Condition,levels=c("SHAM","MI"))
png("celltype.proportion.byCTST.png",width=1500,height=1500,res=300)
p1 <- ggplot(celltype.stats,aes(x=Condition,y=percent))+geom_bar(stat="identity",position="stack",aes(fill=CTST),width=0.75)+theme_bw()+facet_grid(.~Age,space="free",scale="free")
p1 <- p1 + theme(panel.grid=element_blank(),legend.position="bottom")
p1 <- p1 + xlab("") + ylab("Percentage")
print(p1)
dev.off()


write.table(celltype.stats,file="celltype.stats.txt",sep="\t",quote=FALSE,row.names=F)




## extract norm count
## cluster the genes
## divide them and plot using GGPLOT2

##### Back up













### Preparation for cluster identification
metadata = heart.combined@meta.data
metadata$cell =rownames(metadata)

scale.data = GetAssayData(heart.combined,slot = "scale.data")
scaled.data= as.matrix(scale.data)
scaled.data = (scaled.data-rowMeans(scaled.data))/apply(scaled.data,1,sd)
head(scaled.data[,1:20])

plotdata = melt(as.matrix(scaled.data))
names(plotdata)=c("gene","sample","z")
plotdata$cluster = Idents(heart.combined)[match(plotdata$sample,names(Idents(heart.combined)))]
plotdata$condition = metadata$condition[match(plotdata$sample,metadata$cell)]

png("MARKER_featureplot.png",width=6000,height=6000,res=300)
FeaturePlot(heart.combined, features=c("Myl2","Tnnt2","Myh6","Nppa","Dach1","Emcn","Egfl7","Vwf","Cdh5","Tie1","Postn","Dcn","Fn1","Rgs5","Abcc9","Pcdh7","Cd163","Mrc1","Ikzf1","Fyb1"),raster=FALSE)
dev.off()


### Count proportion of cell type
require(stringr)
require(ggplot2)
options(bitmapType="cairo")
celltype.stats = as.data.frame(table(Idents(heart)))
celltype.stats = cbind(celltype.stats,as.data.frame(matrix(unlist(str_split(celltype.stats$Var1,"_")),byrow=TRUE,ncol=3)))
names(celltype.stats)=c("ID","Freq","CT","Condition","Age")
celltype.stats$id = paste(celltype.stats$Condition,celltype.stats$Age)

total = aggregate(celltype.stats$Freq,by=list(celltype.stats$Condition,celltype.stats$Age),FUN=sum)
names(total) = c("Condition","Age","Total")
total$id = paste(total$Condition,total$Age)

celltype.stats$Total = total$Total[match(celltype.stats$id,total$id)]
celltype.stats$percent = celltype.stats$Freq/celltype.stats$Total*100

png("celltype.proportion.byCT.png",width=1500,height=1000,res=300)
p1 <- ggplot(celltype.stats,aes(x=Condition,y=percent))+geom_bar(stat="identity",position="stack",aes(fill=CT),width=0.75)+theme_bw()
p1 <- p1 + theme(panel.grid=element_blank(),legend.position="bottom") + scale_fill_manual(values=c("orangered3","royalblue"))
p1 <- p1 + xlab("") + ylab("Percentage") + coord_flip()
print(p1)
dev.off()



total.bytype = aggregate(celltype.stats$Freq,by=list(celltype.stats$CT),FUN=sum)
names(total.bytype) = c("CT","Total")
celltype.stats$Total.byCT = total.bytype$Total[match(celltype.stats$CT,total.bytype$CT)]
celltype.stats$percent.byCT = celltype.stats$Freq/celltype.stats$Total.byCT*100
write.table(celltype.stats,file="celltype.proportion.txt",sep="\t",quote=FALSE,row.names=F)

temp = celltype.stats[ celltype.stats$Condition == "MI", ]
celltype.stats$rank = temp$percent.byCT[match(celltype.stats$CT,temp$CT)]
temp

png("celltype.proportion.byCT.png",width=1500,height=1000,res=300)
p1 <- ggplot(celltype.stats,aes(x=reorder(CT,rank),y=percent.byCT))+geom_bar(stat="identity",position="stack",aes(fill=Condition),width=0.75)+theme_bw()
p1 <- p1 + theme(panel.grid=element_blank(),legend.position="bottom") + scale_fill_manual(values=c("orangered3","royalblue"))
p1 <- p1 + xlab("") + ylab("Percentage") + coord_flip()
print(p1)
dev.off()





meta = as.data.frame(heart.combined@meta.data)
umap.coord = as.data.frame(heart.combined@reductions$umap@cell.embeddings)
umap.coord$group = meta$condition
umap.coord$group = factor(umap.coord$group,levels=c("SHAM_P1","MI_P1"))
png("UMAP_SHAM_MI_bycondition_P1.png",width=3000,height=3000,res=300)
p1 <- ggplot(umap.coord,aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color=group),size=1.2)+theme_bw()+theme(panel.grid=element_blank())
p1 <- p1 + scale_color_manual(values=c("blue","orangered3"))
print(p1)
dev.off()




