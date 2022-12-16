source('R/lib.R')
md <- readRDS(file = sprintf('results/all_samples_maturation_trajectory_meta_data.Rds'))
cm <- readRDS('data/dropseq_digitial_expression.Rds')

md <- subset(md, postmitotic & eminence == 'MGE')
cm <- cm[, rownames(md)]
cat('There are', ncol(cm), 'post-mitotic cells\n')
pcg <- read.table('annotation/Mus_musculus.GRCm38.84.protein_coding_genes.txt', stringsAsFactors=FALSE)$V1
cm <- cm[rownames(cm) %in% pcg, ]
genes <- rownames(cm)[apply(cm>0, 1, sum) >= 5]
cat('Normalizing', length(genes), 'genes that have been detected in at least 5 post-mitotic cells\n')

if (length(unique(md$sample.name)) > 1) {
  expr <- norm.nb.reg(cm[genes, ], md[, c('reads', 'mols.per.gene', 'cc', 'sample.name')], min.theta=0.01, pr.th=30, bins=64)
} else {
  expr <- norm.nb.reg(cm[genes, ], md[, c('reads', 'mols.per.gene', 'cc')], min.theta=0.01, pr.th=30, bins=64)
}
vg <- genes[scale(sqrt(apply(expr^2, 1, sum)))[, 1] > 1]

vg.expr <- expr[vg, ]
dim(vg.expr)

sample_metadata <- as.data.frame(stringr::str_split_fixed(colnames(vg.expr), "_", 2))
colnames(sample_metadata) <- c("Cell")
sample_metadata$SampleID <- colnames(vg.expr)
colData <- DataFrame(sample_metadata)

gene_metadata <- data.frame(geneName = rownames(vg.expr))
rownames(gene_metadata) <- rownames(vg.expr)
gene_metadata <- DataFrame(gene_metadata)

se <- SummarizedExperiment(
  assays = SimpleList(fpkm = as.matrix(vg.expr)),
  rowData = gene_metadata,
  colData = sample_metadata
)

rowData(se)$geneName
colData(se)$Cell

assay(se,'logfpkm')<-log10(assay(se,'fpkm')+1)
clustData<-assay(se,'fpkm')
nClusters<-3
clustData.scaled<-t(scale(t(clustData)))

clustData.kmeans<-kmeans(clustData.scaled,nClusters)
clustData.plot<-clustData.scaled
clustData.plot<-as.data.frame(clustData.plot)
clustData.plot$gene<-rownames(clustData.plot)
clustData.plot$cluster<-clustData.kmeans$cluster
clustData.melt<-reshape2::melt(clustData.plot,id.vars=c("gene","cluster"))
clustData.melt<-merge(clustData.melt,colData(se),by.x="variable",by.y="SampleID")
clustData.melt<-as.data.frame(clustData.melt)
kmeans.summary <- clustData.melt %>%
  group_by(cluster,Cell) %>%
  summarise(mean=mean(value),sd=sd(value))
#p<-ggplot(kmeans.summary,aes(x=Cell, y=mean, ymax=mean+sd, ymin=mean-sd)) + geom_point(aes(color=Cell)) + geom_line(aes(group=Cell,color=Cell)) + geom_hline(yintercept=0,linetype="dashed") + theme_minimal()
#p + facet_wrap('cluster',scales="free_y")

library(uwot)

genes.umap<-uwot::umap(clustData.scaled,verbose=T)
plot(genes.umap, col=clustData.kmeans$cluster, pch=20)
