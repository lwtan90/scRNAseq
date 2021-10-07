require(Seurat)
library(ggplot2)
library(cowplot)
require(reshape)

#### Subcluster all CM changes for P9
readSIG <- function(filename)
{
	data = read.table(filename,header=TRUE)
	data = data[ data$p_val_adj<0.05, ]
	data$gene = rownames(data)
	data$class = gsub("__P9_DE.txt","",rep(filename))
	return(data)
}

CM0 = readSIG("CM_0__P9_DE.txt")
CM15 = readSIG("CM_15__P9_DE.txt")
CM19 = readSIG("CM_19__P9_DE.txt")
CM1 = readSIG("CM_1__P9_DE.txt")
CM2 = readSIG("CM_2__P9_DE.txt")
CM5 = readSIG("CM_5__P9_DE.txt")
CM7 = readSIG("CM_7__P9_DE.txt")
dCM8 = readSIG("dCM_8__P9_DE.txt")


## hierarchical cluster the cells based on DE
merged = rbind(CM0,CM1,CM2,CM5,CM7,dCM8,CM15,CM19)
as.data.frame(table(merged$gene))
statistics = as.data.frame(table(merged$gene))

heart = readRDS("defined.celltype.RDS")
DefaultAssay(heart) <- "RNA"

info = as.data.frame((matrix(unlist(str_split(heart$condition,"_")),byrow=TRUE,ncol=2)))
names(info)=c("treatment","age") 
heart@meta.data$age = info$age


## subset CM
CM = subset(x = heart, idents = c("CM","dCM"))
CM = subset(x = CM, subset = age == "P9")

scale.data = as.data.frame(CM@assays$RNA@counts)
target.scaledata = scale.data[ rownames(scale.data) %in% merged$gene, ]


### Form average Expression
heart = readRDS("../defined.celltype.RDS")
heart$subtype = paste(Idents(heart),heart$seurat_clusters,heart$condition,sep="_")
heart$CT = Idents(heart)
Idents(heart) <- "subtype"

## calculate gene expression
avg.heart.cells <- as.data.frame(log1p(AverageExpression(heart, verbose = FALSE)$RNA))
r = cor(avg.heart.cells,method="pearson")
options(bitmapType="cairo")
png("correlation.png",width=3500,height=3500,res=300)
heatmap.2(as.matrix(r),Rowv=as.dendrogram(hclust(as.dist(1-r),method="ward.D2")), Colv=as.dendrogram(hclust(as.dist(1-r),method="ward.D2")),col=redgreen(100),trace="none",scale="none",density.info="none")
dev.off()

## just focus on dCM and CM
avg.CM.cells = avg.heart.cells[,grep("CM",colnames(avg.heart.cells))]
r = cor(avg.CM.cells,method="pearson")
options(bitmapType="cairo")
png("correlation_CM.png",width=3500,height=3500,res=300)
heatmap.2(as.matrix(r),Rowv=as.dendrogram(hclust(as.dist(1-r),method="ward.D2")), Colv=as.dendrogram(hclust(as.dist(1-r),method="ward.D2")),col=redgreen(100),trace="none",scale="none",density.info="none")
dev.off()


## Rope in DE genes
readSIG <- function(filename)
{
	data = read.table(filename,header=TRUE)
	data = data[ data$p_val_adj<0.01, ]
	data$gene = rownames(data)
	data$class = gsub("__P9_DE.txt","",rep(filename))
	return(data)
}

CM0 = readSIG("../CM_0__P9_DE.txt")
CM15 = readSIG("../CM_15__P9_DE.txt")
##CM19 = readSIG("../CM_19__P9_DE.txt")
CM1 = readSIG("../CM_1__P9_DE.txt")
CM2 = readSIG("../CM_2__P9_DE.txt")
CM5 = readSIG("../CM_5__P9_DE.txt")
CM7 = readSIG("../CM_7__P9_DE.txt")
dCM8 = readSIG("../dCM_8__P9_DE.txt")


## hierarchical cluster the cells based on DE
merged = rbind(CM0,CM1,CM2,CM5,CM7,dCM8,CM15)
statistics = as.data.frame(table(merged$gene))

DE.CM = avg.CM.cells[ rownames(avg.CM.cells) %in% merged$gene, ]
r = cor(DE.CM,method="pearson")
options(bitmapType="cairo")
png("correlation_CM_DE.png",width=3500,height=3500,res=300)
heatmap.2(as.matrix(r),Rowv=as.dendrogram(hclust(as.dist(1-r),method="ward.D2")), Colv=as.dendrogram(hclust(as.dist(1-r),method="ward.D2")),col=redgreen(100),trace="none",scale="none",density.info="none")
dev.off()

## cluster
hc = hclust(as.dist(1-cor(DE.CM,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(DE.CM),method="spearman")),method="ward.D2")
z = (DE.CM-rowMeans(DE.CM))/apply(DE.CM,1,sd)
z[ z>2] = 2
z[ z<(-2)] = -2
png("heatmap_CM_DE.png",width=3500,height=3500,res=300)
heatmap.2(as.matrix(z),Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),col=redgreen(100),trace="none",scale="row",density.info="none")
dev.off()

## find genes
genetree = as.data.frame(cutree(hr,k=10))
names(genetree)="k"
colors = c("red","orange","yellow","green","blue","magenta","black","grey60","turquoise","brown")
genetree$color = colors[genetree$k]

png("heatmap_CM_DE.png",width=3500,height=3500,res=300)
heatmap.2(as.matrix(z),Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),col=redgreen(100),trace="none",scale="row",density.info="none",RowSideColors=as.character(genetree$color))
dev.off()


## Extract CM-specific single-cell profiles
CM <- subset(heart, subset = CT == c("CM","dCM"))
CM.cells = as.data.frame(heart@assays$RNA@data)
DE.cells = CM.cells[ rownames(CM.cells) %in% merged$gene, ]
z = (DE.cells - rowMeans(DE.cells))/apply(DE.cells,1,sd)
plotdata = melt(as.matrix(z))
names(plotdata)=c("gene","sample","z")
plotdata$z[ plotdata$z>2 ] = 2
plotdata$z[ plotdata$z< (-2) ] = -2
plotdata$genegroup = genetree$color[match(plotdata$gene,rownames(genetree))]
plotdata$group = Idents(CM)[match(plotdata$sample,names(Idents(CM)))]

"dCM_8_SHAM_P1" "CM_0_SHAM_P1"  "CM_5_SHAM_P1"  "CM_1_SHAM_P1" 
"CM_7_SHAM_P1"  "CM_19_SHAM_P1" "CM_2_SHAM_P1"  "CM_15_SHAM_P1"
"CM_1_MI_P1"    "CM_0_MI_P1"    "dCM_8_MI_P1"   "CM_2_MI_P1"   
"CM_5_MI_P1"    "CM_7_MI_P1"    "CM_19_MI_P1"   "CM_15_MI_P1"  
"CM_1_SHAM_P9"  "CM_15_SHAM_P9" "CM_0_SHAM_P9"  "CM_5_SHAM_P9" 
"dCM_8_SHAM_P9" "CM_2_SHAM_P9"  "CM_7_SHAM_P9"  "CM_19_SHAM_P9"
"CM_0_MI_P9"    "CM_5_MI_P9"    "CM_1_MI_P9"    "CM_2_MI_P9"   
["CM_15_MI_P9"   "dCM_8_MI_P9"   "CM_19_MI_P9"   "CM_7_MI_P9"   

## to speed up
plotdata2 = plotdata[ plotdata$sample %in% sample(colnames(DE.cells),ncol(DE.cells)/4), ]
png("CM_masterheatmap.png",width=3000,height=3000,res=300)
ggplot(plotdata2,aes(x=sample,y=gene))+geom_tile(aes(fill=z))+theme_bw()+theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+scale_fill_gradientn(values=rescale(c(-1,0,1)),colors=c("darkblue","turquoise","white","orange","orangered3")
)+facet_grid(genegroup~group,space="free",scale="free")+theme(panel.spacing=unit(0,"lines"))
dev.off()


## Calculate Average CM profiles
## Purpose: Identify which cell types are more similar
CM1 <- subset(CM, idents = "CD4 Naive T")
Idents(t.cells) <- "stim"
avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p1 + p2