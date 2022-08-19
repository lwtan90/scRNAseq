require(ggplot2)

options(bitmapType="cairo")

## meta file 1
## meta file 2
meta1 = read.table("../control/control/seurat/QC.txt",header=TRUE)
meta1$sample = rep("control")


meta2 = read.table("../LMNA1/control/seurat/QC.txt",header=TRUE)
meta2$sample = rep("LMNA1")

filterDATA <- function(testdata)
{
	print(table(testdata$percent.mt>40))
}
filterDATA(meta1)
filterDATA(meta2)

meta = rbind(meta1,meta2)
summary(meta)

#orig.ident	nCount_RNA	nFeature_RNA	nCount_ATAC	nFeature_ATAC	nucleosome_signal	nucleosome_percentile	TSS.enrichment	TSS.percentile	percent.mt

## plot nCount_RNA
p1 <- ggplot(meta,aes(x=nCount_RNA,fill=sample)) + geom_density(alpha=0.6) + theme_bw() + theme(panel.grid=element_blank())
p1 <- p1 + scale_fill_manual(values=c("orange","turquoise")) + scale_x_log10() + geom_vline(xintercept=c(500,100000))
png("nCount_RNA.png",width=1500,height=1300,res=300)
print(p1)
dev.off()

## plot nFeature_RNA
p1 <- ggplot(meta,aes(x=nFeature_RNA,fill=sample)) + geom_density(alpha=0.6) + theme_bw() + theme(panel.grid=element_blank())
p1 <- p1 + scale_fill_manual(values=c("orange","turquoise")) + scale_x_log10() + geom_vline(xintercept=c(400,10000))
png("nFeature_RNA.png",width=1500,height=1300,res=300)
print(p1)
dev.off() 

## plot nCount_ATAC
p1 <- ggplot(meta,aes(x=nCount_ATAC,fill=sample)) + geom_density(alpha=0.6) + theme_bw() + theme(panel.grid=element_blank())
p1 <- p1 + scale_fill_manual(values=c("orange","turquoise")) + scale_x_log10() + geom_vline(xintercept=c(400,90000))
png("nCount_ATAC.png",width=1500,height=1300,res=300)
print(p1)
dev.off()

## plot nFeature_ATAC
p1 <- ggplot(meta,aes(x=nFeature_ATAC,fill=sample)) + geom_density(alpha=0.6) + theme_bw() + theme(panel.grid=element_blank())
p1 <- p1 + scale_fill_manual(values=c("orange","turquoise")) + scale_x_log10() + geom_vline(xintercept=c(500,100000))
png("nFeature_ATAC.png",width=1500,height=1300,res=300)
print(p1)
dev.off()

## plot nucleosome_signal
p1 <- ggplot(meta,aes(x=nucleosome_signal,fill=sample)) + geom_density(alpha=0.6) + theme_bw() + theme(panel.grid=element_blank())
p1 <- p1 + scale_fill_manual(values=c("orange","turquoise"))+geom_vline(xintercept=4)
png("nucleosome_signal.png",width=1500,height=1300,res=300)
print(p1)
dev.off()

## plot TSS.enrichment
p1 <- ggplot(meta,aes(x=TSS.enrichment,fill=sample)) + geom_density(alpha=0.6) + theme_bw() + theme(panel.grid=element_blank())+xlim(0,10)
p1 <- p1 + scale_fill_manual(values=c("orange","turquoise"))+geom_vline(xintercept=2)
png("TSS.enrichment.png",width=1500,height=1300,res=300)
print(p1)
dev.off()

## plot percent.mt
p1 <- ggplot(meta,aes(x=percent.mt,fill=sample)) + geom_density(alpha=0.6) + theme_bw() + theme(panel.grid=element_blank())
p1 <- p1 + scale_fill_manual(values=c("orange","turquoise"))  + geom_vline(xintercept=40)
png("percent.mt.png",width=1500,height=1300,res=300)
print(p1)
dev.off()


## plot correlation between nCount_RNA and nFeature_RNA
p1 <- ggplot(meta,aes(x=nCount_RNA,y=nFeature_RNA)) + geom_bin_2d(bins=1000) +scale_fill_continuous(type = "viridis")+ theme_bw() + theme(panel.grid=element_blank())+scale_x_log10()
p1 <- p1 + facet_grid(.~sample)
png("nCount_RNA_nFeature_RNA.png",width=2800,height=1300,res=300)
print(p1)
dev.off()

## plot correlation between nCount_ATAC and nFeature_ATAC
p1 <- ggplot(meta,aes(x=nCount_ATAC,y=nFeature_ATAC)) + geom_bin_2d(bins=1000) +scale_fill_continuous(type = "viridis")+ theme_bw() + theme(panel.grid=element_blank()) + scale_y_log10()
p1 <- p1 + facet_grid(.~sample)
png("nCount_ATAC_nFeature_ATAC.png",width=2800,height=1300,res=300)
print(p1)
dev.off()

