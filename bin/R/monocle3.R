### This script allows you to run Monocle3 using the raw count from Seurat
## Please edit the following steps manually:
## 1. P1
## 2. P2 (optional)
## 3. P3
### You will need to have monocle3 and Seurat installed
require(Seurat)
require(monocle3)

options(bitmapType="cairo")

# P1
## Edit the R data. If you use RData (rather than RDS), just change to load("CM.RData"), and pass the variable to "heart".
### read in seurat object
heart = readRDS("CM.RDS")
heart$group = Idents(heart)

### Subset cluster 1,2,3,7,8
##heart = subset(heart, subset = group %in% c(2,7,8))
DefaultAssay(heart) = "SCT"

# P2
# If you need to subset the cell type / subgroup, use the following lines.
# You might need to rescale (optional)
### Subset cluster 1,2,3,7,8
##heart = subset(heart, subset = group %in% c(2,7,8))


## read in filtered datasets from Seurat object
## No editting required
# Data refers to the count matrix (here I used normalized count)
data <- as(as.matrix(heart@assays$SCT@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = heart@meta.data)
# Meta.data refers to the metadata requied to form singlecell object for monocle3
meta.data = heart@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## clearing memory
rm(heart)

## run monocle
monocle_cds <- new_cell_data_set(data,cell_metadata=meta.data,gene_metadata=fData)
monocle_cds <- preprocess_cds(monocle_cds, num_dim = 100)
saveRDS(monocle_cds,file="preprocess.cds.RDS")
monocle_cds <- reduce_dimension(monocle_cds,reduction_method="UMAP")
monocle_cds <- cluster_cells(monocle_cds)
monocle_cds <- learn_graph(monocle_cds,use_partition = F)
monocle_cds$subtype = monocle_cds$cluster
saveRDS(monocle_cds,file="graph.learned.nopartition.RDS")

# P3
## Depending on how you run, if you use RStudio, there might be an interface for you to pick the root.
## Else, please edit the "time_bin=4" to the specific subtype that you think is the starting point of the trajectory
## Next, change the "group" to the corresponding cell type column in the metadata. Here i used group because that was how I set up my metadata. It can be changed to "celltype" or "cluster" or etc.
## select root
cds = monocle_cds
get_earliest_principal_node <- function(cds, time_bin="4"){
  cell_ids <- which(colData(monocle_cds)[, "group"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

## congratulation. you have your pseudotime
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
saveRDS(cds,file="trajectory.100.nopartition.RDS")

## plot trajectory
png("tracjetory.nopartition.png",width=2500,height=2500,res=300)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
dev.off()

png("tracjetory.nopartition.labelBRANCH.png",width=2500,height=2500,res=300)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=TRUE,graph_label_size=1.5)
dev.off()

## optional
cds$group = paste("G",cds$group)
png("tracjetory.nopartition.coloredsubgroup.png",width=2500,height=2500,res=300)
plot_cells(cds,color_cells_by = "group",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=TRUE,graph_label_size=1.5)
dev.off()
	   

