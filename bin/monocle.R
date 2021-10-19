### Merging seurat with Monocle
require(Seurat)
require(monocle3)

#Load Seurat object
seurat_object <- readRDS('CM.RDS')

## Extract degenes only
##degenes = read.table("geneclass.txt",header=TRUE)
##CM = subset(seurat_object,features = rownames(degenes))
degenes = read.table("CM.DE.txt")
CM = subset(seurat_object,features = degenes$V1)

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(CM@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = CM@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- new_cell_data_set(data,cell_metadata=CM@meta.data,gene_metadata=fData)
monocle_cds <- preprocess_cds(monocle_cds, num_dim = 100)

## Step 3: Reduce the dimensions using UMAP
monocle_cds <- reduce_dimension(monocle_cds)

## Step 4: Cluster the cells
monocle_cds <- cluster_cells(monocle_cds)

## Step 5: Learn a graph
monocle_cds <- learn_graph(monocle_cds)

monocle_cds$group = paste(monocle_cds$condition,monocle_cds$integrated_snn_res.1,sep="_")

### subtype1
# a helper function to identify the root principal points:
cds = monocle_cds
get_earliest_principal_node <- function(cds, time_bin="SHAM_P1_0"){
  cell_ids <- which(colData(monocle_cds)[, "group"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

options(bitmapType="cairo")
png("tracjetory_CM.png",width=2500,height=2500,res=300)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
dev.off()

png("tracjetory_CMsubtype.png",width=2500,height=2500,res=300)
plot_cells(cds,color_cells_by = "group",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
dev.off()

png("tracjetory_condition.png",width=2500,height=2500,res=300)
plot_cells(cds,color_cells_by = "condition",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
dev.off()
CM$pseudotime = cds$pseudotime
saveRDS(CM,file="pseudoadded_CM.RDS")

## incase
original.cds = cds

cds$pseudotime = pseudotime(cds)
cds_subset = cds[ cds$pseudotime!=Inf, ]

genegroup = read.table("geneclass.txt",header=TRUE)
cds_subset <- cds_subset[rowData(cds_subset)$gene_short_name %in% rownames(genegroup),]

### Unfortunaley 
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
cds_subset_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
write.table(cds_subset_test_res,file="pseudotime_graphtest.txt",sep="\t",quote=FALSE)

options(bitmapType="cairo")
png("genemarker_module.png",width=3000,height=2500,res=300)
plot_cells(cds_subset, genes=head(rownames(cds_subset_test_res),12),show_trajectory_graph=FALSE,label_cell_groups=FALSE, label_leaves=FALSE)
dev.off()

		   
cds_subset_deg_ids <- row.names(subset(cds_subset_test_res, q_value < 0.01))
### Generating gene modules
## try increating resoltion
gene_module_df <- find_gene_modules(cds_subset, resolution=1e-1)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)), cell_group=cds_subset$subtype)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

png("module_CM.png",width=2000,height=2000,res=300)
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,scale="column", clustering_method="ward.D2",fontsize=6)
dev.off()

png("plotcell_CM.png",width=3000,height=3000,res=300)
plot_cells(cds_subset, genes=gene_module_df %>% filter(module %in% c(1, 2, 3, 4,5,6)),group_cells_by="cluster",color_cells_by="cluster",show_trajectory_graph=FALSE)
dev.off()





## Plotting the outcome
## Plotting the outcome
CM = readRDS("pseudoadded_CM.RDS")
CM = subset(CM, subset = pseudotime !=Inf)
genegroup = read.table("gene_module_left.txt",header=TRUE)

#Set 1: Pseudotime of the cardiomyocytes
require(ggplot2)
options(bitmapType="cairo")
metadata = CM@meta.data

metadata$treatment = gsub("_P[0-9]","",metadata$condition)
metadata$time = gsub("SHAM_","",metadata$condition)
metadata$time = gsub("MI_","",metadata$time)


png("pseudotime_celltype_left.png",width=2000,height=3500,res=300)
p1<-ggplot(metadata,aes(x=reorder(subtype,pseudotime),y=pseudotime))+geom_boxplot(aes(fill=time))+theme_bw()+theme(panel.grid=element_blank())
p1<-p1 + scale_fill_manual(values=c("green","blue"))+ coord_flip()
p1 <- p1 + facet_grid(treatment~.,space="free",scale="free")
print(p1)
dev.off()

#Set 2: Average module expression vs pseudotime
### extracting cell counts
CM.cells = as.data.frame(CM@assays$RNA@data)
DE.cells = CM.cells[ rownames(CM.cells) %in% genegroup$id, ]

## z-norm
z = (DE.cells - rowMeans(DE.cells))/apply(DE.cells,1,sd)

### form plotdata
require(reshape)
plotdata = melt(as.matrix(z))
names(plotdata)=c("gene","sample","z")
plotdata$genegroup = genegroup$module[match(plotdata$gene,genegroup$id)]
plotdata$group = CM$pseudotime[match(plotdata$sample,names(Idents(CM)))]

genevariable = aggregate(plotdata$z,by=list(plotdata$sample,plotdata$genegroup),FUN=mean)
names(genevariable)=c("cell","genegroup","z")

##library(reshape2)
#genevariable.df = acast(genevariable, cell~genegroup, value.var="z")
#hc = hclust(as.dist(1-cor(genevariable.df)),method="ward.D2")
#require(gplots)
#png("heatmap_modulegenegroup.png",width=1500,height=2500,res=300)
##heatmap.2(as.matrix(genevariable.df),Rowv=FALSE,Colv=as.dendrogram(hc),scale="none",trace="none",density.info="none",col=redgreen(100),labRow="none")
#dev.off()


genevariable$pseudotime = CM$pseudotime[match(genevariable$cell,names(Idents(CM)))]

cor.df = data.frame(genegroup=c(),r=c())
for(i in 1:max(genevariable$genegroup))
{
	x = genevariable$z[ genevariable$genegroup == i ]
	y = genevariable$pseudotime[ genevariable$genegroup == i ]
	print(cor(x,y))
	cor.df = rbind(cor.df,data.frame(genegroup=c(i),r=cor(x,y)))
}

genevariable$rank = cor.df$r[match(genevariable$genegroup,cor.df$genegroup)]
##plotdata$z[ plotdata$z>2 ] = 2
##plotdata$z[ plotdata$z< (-2) ] = -2

require(ggplot2)
require(scales)
options(bitmapType="cairo")
genevariable2 = genevariable
genevariable2$z[ genevariable2$z>0.5]=0.5
genevariable2$z[ genevariable2$z<(-0.5)]=-0.5
png("pseudotime_CM_left.png",width=3000,height=1500,res=300)
p1<-ggplot(genevariable2,aes(x=reorder(cell,pseudotime),y=reorder(as.factor(genegroup),rank)))+geom_tile(aes(fill=z))+theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank())
p1<-p1+scale_fill_gradientn(values=rescale(c(-0.5,0,0.5)),colors=c("blue","turquoise","white","orange","red","darkred"))+theme(panel.spacing=unit(0,"lines"))
print(p1)
dev.off()

genevariable$subtype = metadata$subtype[match(genevariable$cell,rownames(metadata))]
genevariable$condition = gsub("CM_[0-9]_","", genevariable$subtype)
genevariable$condition = gsub("CM_[0-9][0-9]_","", genevariable$condition)
genevariable$condition = gsub("d","", genevariable$condition)
genevariable$ST = gsub("CM_","",genevariable$subtype)
genevariable$ST = gsub("_SHAM_P[0-9]","",genevariable$ST)
genevariable$ST = gsub("_MI_P[0-9]","",genevariable$ST)
genevariable$condition = factor(genevariable$condition,levels=c("SHAM_P1","MI_P1","SHAM_P9","MI_P9"))
genevariable2 = genevariable
genevariable2$z[ genevariable2$z>0.5]=0.5
genevariable2$z[ genevariable2$z<(-0.5)]=-0.5
png("pseudotime_CM_left_condition.png",width=3000,height=1500,res=300)
p1<-ggplot(genevariable2,aes(x=reorder(cell,pseudotime),y=reorder(as.factor(genegroup),rank)))+geom_tile(aes(fill=z))+theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank())
p1<-p1+scale_fill_gradientn(values=rescale(c(-0.5,0,0.5)),colors=c("blue","turquoise","white","orange","red","darkred"))+theme(panel.spacing=unit(0,"lines"))
p1<-p1+facet_grid(.~condition,space="free_y",scale="free")
print(p1)
dev.off()

### plot module-specific expression
genevariable$subtype = metadata$subtype[match(genevariable$cell,rownames(metadata))]
genevariable$condition = gsub("CM_[0-9]_","", genevariable$subtype)
genevariable$condition = gsub("CM_[0-9][0-9]_","", genevariable$condition)
genevariable$condition = gsub("d","", genevariable$condition)
genevariable$ST = gsub("CM_","",genevariable$subtype)
genevariable$ST = gsub("_SHAM_P[0-9]","",genevariable$ST)
genevariable$ST = gsub("_MI_P[0-9]","",genevariable$ST)
genevariable$condition = factor(genevariable$condition,levels=c("SHAM_P1","MI_P1","SHAM_P9","MI_P9"))
png("expression_subtype_module_left.png",width=3000,height=3500,res=300)
p1<-ggplot(genevariable,aes(x=as.factor(condition),y=z))+geom_violin(trim=FALSE,aes(fill=condition))+theme_bw()+theme(panel.grid=element_blank())
p1<-p1+theme(panel.spacing=unit(0,"lines"))+scale_fill_manual(values=c("lightblue","orange","blue","red"))
p1<-p1+facet_grid(genegroup~ST,space="free_x",scale="free")
print(p1)
dev.off()

for(i in 1:max(genevariable$genegroup))
{
	genevariable3 = genevariable[ genevariable$genegroup==i , ]
	png(paste("expression_left_subtype",i,".png",sep=""),width=3000,height=3500,res=300)
	p1<-ggplot(genevariable3,aes(x=pseudotime,y=z))+geom_point(aes(color=condition))+theme_bw()+theme(panel.grid=element_blank())
	p1<-p1+theme(panel.spacing=unit(0,"lines"))+scale_color_manual(values=c("lightblue","orange","blue","red"))
	p1<-p1+facet_grid(ST~.,space="free",scale="free_x")
	p1<-p1+stat_smooth(method="loess",se=FALSE,span=0.75)
	print(p1)
	dev.off()
}

#####







## right
cds = monocle_cds
get_earliest_principal_node <- function(cds, time_bin="CM_1_SHAM_P9"){
  cell_ids <- which(colData(monocle_cds)[, "subtype"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

options(bitmapType="cairo")
png("tracjetory_CM_right.png",width=2500,height=2500,res=300)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
dev.off()

png("tracjetory_CMsubtype_right.png",width=2500,height=2500,res=300)
plot_cells(cds,color_cells_by = "subtype",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
dev.off()

CM$pseudotime2 = pseudotime(cds)
saveRDS(CM,file="pseudoadded_CM_right.RDS")

## incase
original.cds = cds

cds$pseudotime = pseudotime(cds)
cds_subset = cds

### Unfortunaley 
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
cds_subset_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
write.table(cds_subset_test_res,file="pseudotime_graphtest_right.txt",sep="\t",quote=FALSE)

options(bitmapType="cairo")
png("genemarker_module_right.png",width=3000,height=2500,res=300)
plot_cells(cds_subset, genes=head(rownames(cds_subset_test_res),12),show_trajectory_graph=FALSE,label_cell_groups=FALSE, label_leaves=FALSE)
dev.off()

		   
cds_subset_deg_ids <- row.names(subset(cds_subset_test_res, q_value < 0.01))
### Generating gene modules
## try increating resoltion
gene_module_df <- find_gene_modules(cds_subset, resolution=1e-1)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)), cell_group=cds_subset$subtype)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

png("module_CM_right.png",width=2000,height=2000,res=300)
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,scale="column", clustering_method="ward.D2",fontsize=6)
dev.off()

png("plotcell_CM_right.png",width=3000,height=3000,res=300)
plot_cells(cds_subset, genes=gene_module_df %>% filter(module %in% c(1, 2, 3, 4,5,6,7,8,9,10)),group_cells_by="cluster",color_cells_by="cluster",show_trajectory_graph=FALSE)
dev.off()

write.table(gene_module_df,file="gene_module_right.txt",sep="\t",quote=FALSE)

#####



## Plotting the outcome
CM = readRDS("pseudoadded_CM_right.RDS")
CM = subset(CM, subset = pseudotime !=Inf)
genegroup = read.table("gene_module_right.txt",header=TRUE)

#Set 1: Pseudotime of the cardiomyocytes
require(ggplot2)
options(bitmapType="cairo")
metadata = CM@meta.data

metadata$treatment = gsub("_P[0-9]","",metadata$condition)
metadata$time = gsub("SHAM_","",metadata$condition)
metadata$time = gsub("MI_","",metadata$time)


png("pseudotime_celltype_right.png",width=2000,height=3500,res=300)
p1<-ggplot(metadata,aes(x=reorder(subtype,pseudotime),y=pseudotime))+geom_boxplot(aes(fill=time))+theme_bw()+theme(panel.grid=element_blank())
p1<-p1 + scale_fill_manual(values=c("green","blue"))+ coord_flip()
p1 <- p1 + facet_grid(treatment~.,space="free",scale="free")
print(p1)
dev.off()

#Set 2: Average module expression vs pseudotime
### extracting cell counts
CM.cells = as.data.frame(CM@assays$RNA@data)
DE.cells = CM.cells[ rownames(CM.cells) %in% genegroup$id, ]

## z-norm
z = (DE.cells - rowMeans(DE.cells))/apply(DE.cells,1,sd)

### form plotdata
require(reshape)
plotdata = melt(as.matrix(z))
names(plotdata)=c("gene","sample","z")
plotdata$genegroup = genegroup$module[match(plotdata$gene,genegroup$id)]
plotdata$group = CM$pseudotime[match(plotdata$sample,names(Idents(CM)))]

genevariable = aggregate(plotdata$z,by=list(plotdata$sample,plotdata$genegroup),FUN=mean)
names(genevariable)=c("cell","genegroup","z")

##library(reshape2)
#genevariable.df = acast(genevariable, cell~genegroup, value.var="z")
#hc = hclust(as.dist(1-cor(genevariable.df)),method="ward.D2")
#require(gplots)
#png("heatmap_modulegenegroup.png",width=1500,height=2500,res=300)
##heatmap.2(as.matrix(genevariable.df),Rowv=FALSE,Colv=as.dendrogram(hc),scale="none",trace="none",density.info="none",col=redgreen(100),labRow="none")
#dev.off()


genevariable$pseudotime = CM$pseudotime[match(genevariable$cell,names(Idents(CM)))]

cor.df = data.frame(genegroup=c(),r=c())
for(i in 1:max(genevariable$genegroup))
{
	x = genevariable$z[ genevariable$genegroup == i ]
	y = genevariable$pseudotime[ genevariable$genegroup == i ]
	print(cor(x,y))
	cor.df = rbind(cor.df,data.frame(genegroup=c(i),r=cor(x,y)))
}

genevariable$rank = cor.df$r[match(genevariable$genegroup,cor.df$genegroup)]
##plotdata$z[ plotdata$z>2 ] = 2
##plotdata$z[ plotdata$z< (-2) ] = -2

require(ggplot2)
require(scales)
options(bitmapType="cairo")
genevariable2 = genevariable
genevariable2$z[ genevariable2$z>0.5]=0.5
genevariable2$z[ genevariable2$z<(-0.5)]=-0.5
png("pseudotime_CM_right.png",width=3000,height=1500,res=300)
p1<-ggplot(genevariable2,aes(x=reorder(cell,pseudotime),y=reorder(as.factor(genegroup),rank)))+geom_tile(aes(fill=z))+theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank())
p1<-p1+scale_fill_gradientn(values=rescale(c(-0.5,0,0.5)),colors=c("blue","turquoise","white","orange","red","darkred"))+theme(panel.spacing=unit(0,"lines"))
print(p1)
dev.off()

genevariable$subtype = metadata$subtype[match(genevariable$cell,rownames(metadata))]
genevariable$condition = gsub("CM_[0-9]_","", genevariable$subtype)
genevariable$condition = gsub("CM_[0-9][0-9]_","", genevariable$condition)
genevariable$condition = gsub("d","", genevariable$condition)
genevariable$ST = gsub("CM_","",genevariable$subtype)
genevariable$ST = gsub("_SHAM_P[0-9]","",genevariable$ST)
genevariable$ST = gsub("_MI_P[0-9]","",genevariable$ST)
genevariable$condition = factor(genevariable$condition,levels=c("SHAM_P1","MI_P1","SHAM_P9","MI_P9"))
genevariable2 = genevariable
genevariable2$z[ genevariable2$z>0.5]=0.5
genevariable2$z[ genevariable2$z<(-0.5)]=-0.5
png("pseudotime_CM_right_condition.png",width=3000,height=1500,res=300)
p1<-ggplot(genevariable2,aes(x=reorder(cell,pseudotime),y=reorder(as.factor(genegroup),rank)))+geom_tile(aes(fill=z))+theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank())
p1<-p1+scale_fill_gradientn(values=rescale(c(-0.5,0,0.5)),colors=c("blue","turquoise","white","orange","red","darkred"))+theme(panel.spacing=unit(0,"lines"))
p1<-p1+facet_grid(.~condition,space="free_y",scale="free")
print(p1)
dev.off()

### plot module-specific expression
genevariable$subtype = metadata$subtype[match(genevariable$cell,rownames(metadata))]
genevariable$condition = gsub("CM_[0-9]_","", genevariable$subtype)
genevariable$condition = gsub("CM_[0-9][0-9]_","", genevariable$condition)
genevariable$condition = gsub("d","", genevariable$condition)
genevariable$ST = gsub("CM_","",genevariable$subtype)
genevariable$ST = gsub("_SHAM_P[0-9]","",genevariable$ST)
genevariable$ST = gsub("_MI_P[0-9]","",genevariable$ST)
genevariable$condition = factor(genevariable$condition,levels=c("SHAM_P1","MI_P1","SHAM_P9","MI_P9"))
png("expression_subtype_module_right.png",width=3000,height=3500,res=300)
p1<-ggplot(genevariable,aes(x=as.factor(condition),y=z))+geom_violin(trim=FALSE,aes(fill=condition))+theme_bw()+theme(panel.grid=element_blank())
p1<-p1+theme(panel.spacing=unit(0,"lines"))+scale_fill_manual(values=c("lightblue","orange","blue","red"))
p1<-p1+facet_grid(genegroup~ST,space="free_x",scale="free")
print(p1)
dev.off()

for(i in 1:max(genevariable$genegroup))
{
	genevariable3 = genevariable[ genevariable$genegroup==i , ]
	png(paste("expression_right_subtype",i,".png",sep=""),width=3000,height=3500,res=300)
	p1<-ggplot(genevariable3,aes(x=pseudotime,y=z))+geom_point(aes(color=condition))+theme_bw()+theme(panel.grid=element_blank())
	p1<-p1+theme(panel.spacing=unit(0,"lines"))+scale_color_manual(values=c("lightblue","orange","blue","red"))
	p1<-p1+facet_grid(ST~.,space="free",scale="free_x")
	p1<-p1+stat_smooth(method="loess",se=FALSE,span=0.75)
	print(p1)
	dev.off()
}


#####
