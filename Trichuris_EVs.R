library('tximport')
library('DESeq2')
library("pheatmap")
library("vsn")
library("RColorBrewer")
library("ggplot2")
library("reshape2")
library("biomaRt")
library("gProfileR")
library("plyr")
library("gridExtra")
library("EnhancedVolcano")

#paths to data directories
dir <- '~/git_repos/Trichuris_EVs'
tx2gene <- read.table('~/git_repos/Trichuris_EVs/transcripts2genes.E97.tsv', header = T)

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
experiment_design <- read.table(paste0(dir,'/experiment_design.tsv'),header=TRUE)
experiment_design <- experiment_design[which(experiment_design$condition != 'Hpoly_ev'),]
files<-file.path(dir, experiment_design$file,'abundance.h5')
names(files) <- experiment_design$file

#import the table of read numbers and counts and plot
read_counts <- read.table(paste0(dir,'/counts_table.tsv'),header=FALSE)
names(read_counts) <- c("sample","total_reads","kallisto_counts")
read_counts<-melt(read_counts,id.vars = "sample")
p<-ggplot(data=read_counts, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") + 
  theme_minimal() +
  xlab('Sample') + 
  ylab('Count') +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=20), axis.text.x=element_text(angle=90))

pdf(paste0(dir,'/read_counts.pdf'))
print(p)
dev.off()


#import the dataset
txi <- tximport(files, type = "kallisto",tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport (txi,colData=experiment_design,design =~ condition)

#run the analysis
dds <- DESeq(dds)

res <-  results(dds,contrast = c("condition","Trichuris_ev","PBS"),alpha=0.05)
res <- res[order(res$pvalue),]

#export the ordered results files
write.table(res,file=paste0(dir,'/deseq2_results.tsv'),sep="\t",eol="\n")

#biomart function
#for significantly regulated genes, return log fold changes, gene names and gene descriptions.
mart<-function(results_object){
  values <- row.names(subset(results_object, padj < 0.05))
  filters <- ('ensembl_gene_id')
  attributes <- c('ensembl_gene_id','external_gene_name','description')
  query <- getBM(attributes = attributes, filters = filters, values = values, mart = ensembl)
  results_object.df<-as.data.frame(results_object)
  results_object.df<-merge(query,results_object.df,by.x = "ensembl_gene_id",by.y="row.names",all.x=TRUE)
  return(results_object.df)
}

#gprofiler function
#for significantly regulated genes, return 
profile<-function(results_object,sig_results){
  universe <- rownames(subset(results_object, results_object$baseMean > 0))
  genes <- unique(sig_results$ensembl_gene_id)
  gp<-gprofiler(genes, organism="mmusculus", custom_bg = universe)
  return(gp)
}

res_sig<-mart(res)
write.table(res_sig,file=paste0(dir,'/deseq2_sig_genes_mart.tsv'),sep="\t",eol="\n")
res_gp<-profile(res,res_sig)
write.table(res_gp,file = paste0(dir,'/deseq2_sig_genes_gprofiler.tsv'), sep ='\t', eol = '\n')

#Transform the data for visualisation
ntd <- normTransform(dds)
vsd <- vst(dds,blind=FALSE) 
rld <- rlog(dds)

#Choose which transformation to use
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd)) #choose this one
meanSdPlot(assay(rld))


#PCA plot
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
pcaData$condition <- revalue(pcaData$condition, c("Trichuris_ev"="T. muris EVs"))
pcaData$group <- revalue(pcaData$group, c("Trichuris_ev"="T. muris EVs"))
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(paste0(dir,'/PCA.pdf'))
print(ggplot(pcaData, aes(PC1, PC2, color=condition)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove background and grid lines
              text = element_text(size=10)) + 
        coord_fixed() + 
        scale_color_manual(values=c("black","grey")) + 
        labs(color='Treatment')
        ) 
dev.off()

#function to return a counts plot
make_graph<-function(x){
  d <- plotCounts(dds, gene= x[[1]], intgroup="condition", returnData=TRUE)
  d$condition <- revalue(d$condition, c("Hpoly_ev"="H. poly", "Trichuris_ev"="T. muris"))
  d$condition<-factor(d$condition, level = c("PBS","H. poly","T. muris"))
  g<-ggplot(d, aes(x=condition, y=count)) + 
    geom_point(position=position_jitter(w=0.1,h=0),size=2) +
    labs(title = paste0(x[[2]],"-",x[[1]]))
  return(g)
}

#Interesting genes were selected using InnateDB

#make heatmaps of interesting genes
viral_genes<-scan(file = paste0(dir,'/virus.txt'), sep = "\n", what = character())
viral_ids<-res_sig[res_sig$external_gene_name %in% viral_genes, c("ensembl_gene_id","external_gene_name","log2FoldChange" )]
viral_ids<-viral_ids[which(viral_ids$log2FoldChange > )]
values<-as.data.frame(assay(vsd)[viral_ids$ensembl_gene_id,])
row.names(values)<-viral_ids$external_gene_name
names(values)<-c("PBS1","PBS2","PBS3","TmEV1","TmEV2","TmEV3")
colors <- colorRampPalette((brewer.pal(9, "Reds")) )(255)
#pheatmap(values,cluster_cols=F,col=colors)
#Heat map of means
mean_values<-data.frame(PBS = rowMeans(values[,c("PBS1", "PBS2", "PBS3")]), TmEV = rowMeans(values[,c("TmEV1", "TmEV2", "TmEV3")]) )
pdf(paste0(dir,'/viral_heat_map.pdf'))
pheatmap(mean_values,cluster_cols=F,col=colors,treeheight_row = 0, treeheight_col = 0,cellwidth = 20, cellheight = 15)
dev.off()

other_immune<-scan(file = paste0(dir,'/innate_immune_not_viral.txt'), sep = "\n", what = character())
other_immune_ids<-res_sig[res_sig$external_gene_name %in% other_immune, c("ensembl_gene_id","external_gene_name","log2FoldChange" )]
values<-as.data.frame(assay(vsd)[other_immune_ids$ensembl_gene_id,])
row.names(values)<-other_immune_ids$external_gene_name
names(values)<-c("PBS1","PBS2","PBS3","TmEV1","TmEV2","TmEV3")
mean_values<-data.frame(PBS = rowMeans(values[,c("PBS1", "PBS2", "PBS3")]), TmEV = rowMeans(values[,c("TmEV1", "TmEV2", "TmEV3")]) )
pdf(paste0(dir,'/other_immune_heat_map.pdf'))
pheatmap(mean_values,cluster_cols=F,col=colors,treeheight_row = 0, treeheight_col = 0,cellwidth = 20, cellheight = 15)
dev.off()


exosomes<-scan(file = paste0(dir,'/exosome.txt'), sep = "\n", what = character())
exosome_ids<-res_sig[res_sig$external_gene_name %in% exosomes, c("ensembl_gene_id","external_gene_name","log2FoldChange" )]
values<-as.data.frame(assay(vsd)[exosome_ids$ensembl_gene_id,])
row.names(values)<-exosome_ids$external_gene_name
names(values)<-c("PBS1","PBS2","PBS3","TmEV1","TmEV2","TmEV3")
mean_values<-data.frame(PBS = rowMeans(values[,c("PBS1", "PBS2", "PBS3")]), TmEV = rowMeans(values[,c("TmEV1", "TmEV2", "TmEV3")]) )
pdf(paste0(dir,'/exosomes_heat_map.pdf'))
pheatmap(mean_values,cluster_cols=F,col=colors,treeheight_row = 0, treeheight_col = 0,cellwidth = 20, cellheight = 15)
dev.off()


#nice volcano plot
filters <- ('ensembl_gene_id')
attributes <- c('ensembl_gene_id','external_gene_name','description')
query <- getBM(attributes = attributes, filters = filters, values = row.names(res), mart = ensembl)
res.df<-as.data.frame(res)
res.df<-merge(query,res.df,by.x = "ensembl_gene_id",by.y="row.names",all.x=TRUE)
pdf(paste0(dir,'/volcano.pdf'))
EnhancedVolcano(res.df,lab=res.df$external_gene_name,x='log2FoldChange',y='padj',xlim = c(-4,4.5),gridlines.major = FALSE,gridlines.minor = FALSE)
dev.off()

