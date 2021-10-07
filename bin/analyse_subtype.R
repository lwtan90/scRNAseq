
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