require(Seurat)
require(Signac)
require(ggplot2)

options(bitmapType="cairo")

## list directory
top10marker <- function(filename)
{
	data = read.table(filename,header=TRUE)
	data = data[ data$p_val_adj<0.05, ]
	data = data[ order(-data$avg_log2FC),]
	return(head(rownames(data),10))
}

features = c()
file.list = list.files(pattern="^marker_")
for( f in file.list)
{
	print(f)
	features = c(features, top10marker(f))
}

length(unique(features))
features = unique(features)
features = features[ grep("MT-",features,invert=TRUE) ]
features

write.table(features,file="markers.txt",sep="\t",quote=FALSE,row.names=F)

## extract top 5 markers per cluster
## plot the dotplot
heart.combined = readRDS("SCT.prepped.RDS")
##Idents(heart.combined) = heart.combined$integrated_snn_res.0.5
table(Idents(heart.combined))

levels(heart.combined) <- c('CM', 'FB', 'EC', 'MC','Muscle-like','SMC','UNK1','UNK2','UNK3')

p1 <- DotPlot(heart.combined, features = features)
png("top5marker_dotplot.png",height=5000,width=2000,res=300) 
p1 <- p1 & coord_flip() & theme(axis.text.x = element_text(angle = 90))
p1
dev.off()

png("top5marker_heatmap.png",width=5000,height=5000,res=300)
DoHeatmap(heart.combined, features = features, cells = 1:500, size = 4,angle = 90) + NoLegend()
dev.off()

## plot the doheatmap

