library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot);library(ggthemes)
library(DoubletFinder);library(copykat);#library(CaSpER);library(GenomicRanges);library(future)
library(presto);library(harmony);library(ggrepel)
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library(presto);library(harmony);library(ggrepel)
library(clusterProfiler);library(org.Hs.eg.db);library(Rcpp)
library(homologene)
library(KEGGREST)
library(pathview) 
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot);library(ggthemes)
library(DoubletFinder);library(infercnv);library(GenomicRanges);library(future)
library(presto);library(harmony);library(ggrepel)
library(Startrac);library(scRepertoire)
library(clusterProfiler);library(org.Hs.eg.db);library(Rcpp)
library(ROGUE)
library(SingleCellExperiment)
library(smoother)
library(plotrix)
library(ggpubr)
library(RColorBrewer)
library(Seurat)
library(velocyto.R)
library(magrittr)
library(stringr)
library(SeuratDisk)
library(ggalluvial)
library(ggforce)
library(SeuratWrappers)
library(org.Mm.eg.db)
library(AUCell)
library(ape)
library(loomR)
library(Hmisc)
library(irGSEA)
library(CytoTRACE)
library(genomicInstability)
library(Chord)
library(scMetabolism)
library(projectLSI)
library(ggridges)
#BiocManager::install("maftools")
library(maftools)
library(CellChat)
library(ConsensusClusterPlus)
library(loomR)
library(SCopeLoomR)
library(SCENIC)
library(plot1cell)
library(URD)
library(monocle3)
library(tidyr)
library(clusterProfiler)
library(GseaVis)
library(org.Mm.eg.db)

library(enrichplot)
library(tidyverse)
library(ggstatsplot)

set.seed(1234)

main.path='/data4/huanggy/CTTN/'
save.data='/data4/huanggy/CTTN/data/'
save.pic='/data4/huanggy/CTTN/pic/'




#------------CTC_database counts--------loading-----------------------------
#CTC
files = list.files('/data4/huanggy/CTTN/CTC/')
tpm = files[!is.na(str_extract(files,('RSEM.count.+')))]%>%str_remove('RSEM.count.')%>%str_remove('\\.csv')


#Human-CTC-database
human = c("CNP0000095","GSE109761","GSE111065","GSE111842","GSE51827","GSE55807","GSE67939","GSE67980","GSE75367","GSE86978","PRJNA522885","PRJNA603782","PRJNA603789","PRJNA662599")


counts_human = list()
for (i in human) {
  counts_human[[i]] <- read.csv(paste0('/data4/huanggy/CTTN/CTC/RSEM.count.',i,'.csv'),header = T,sep = ',',row.names = 1)
  colnames(counts_human[[i]]) <-  paste0(colnames(counts_human[[i]]),';',i)
}

counts_human = Reduce(cbind,counts_human)
counts_human[,'name'] = (rownames(counts_human) %>% str_split_fixed('_',2))[,2]
counts_human = counts_human[!duplicated(counts_human$name),]
rownames(counts_human) = counts_human$name
counts_human$name = NULL



CTC.human = data.frame(name = colnames(counts_human),source = str_split_fixed(colnames(counts_human),'_',4)[,3],cancertype= str_split_fixed(colnames(counts_human),'_',4)[,1])
rownames(CTC.human) = CTC.human$name
CTC.human[,'patient'] = str_split_fixed(colnames(counts_human),'_',5)[,2]
CTC.human[,'paper'] = str_split_fixed(colnames(counts_human),';',2)[,2]
CTC.human[,'group'] <- paste0(CTC.human[,'source'],'_',CTC.human[,'paper'])
CTC.human[,'group1'] <-  str_split_fixed(colnames(counts_human),'_',5)[,4]

CTC.human[,'group2'] <-  str_split_fixed(str_split_fixed(colnames(counts_human),';',2)[,1],'_',6)[,5]



CTC.human$group <- factor(CTC.human$group,levels = unique(CTC.human$group,ordered = T))



CTC.human1=CTC.human[CTC.human$source %in% c('BloodCTC','Whitecell','WhiteCell'),]
#exclude CTC clusters
CTC.human1=CTC.human1[!CTC.human1$group1%in% c('n'),]
#exclude CTC clusters
CTC.human1=CTC.human1[!CTC.human1$group2%in% c('n'),]
#exclude CTC clusters
CTC.human1=CTC.human1[!CTC.human1$group2%in% c('n'),]





counts_human1 =counts_human[,rownames(CTC.human1[CTC.human1$source %in% c('BloodCTC','Whitecell','WhiteCell'),])]

#-----delete mito & Doublet------

#mito
counts <- CreateSeuratObject(counts =counts_human1)
counts@meta.data$source <-str_split_fixed(rownames(counts@meta.data),'_',4)[,3] 
counts[["percent.mt"]] = PercentageFeatureSet(counts, pattern = "^[Mm]T-")
counts[["RP.mt"]] = PercentageFeatureSet(counts, pattern = "^RP[LS]\\w")
counts <- subset(counts, subset = (nFeature_RNA > 300 & nFeature_RNA< 5000 & percent.mt < 20 & RP.mt < 50))
counts=counts@assays$RNA@counts


#Doublet
  tmp <- CreateSeuratObject(counts) %>% NormalizeData(verbose =F) %>% ScaleData(verbose =F) %>%  FindVariableFeatures(selection.method = "vst", nfeatures = 2000,verbose =F) %>% RunPCA(verbose =F) %>% FindNeighbors(dims =1:10,verbose =F) %>% FindClusters(resolution =0.5,verbose =F)%>% RunUMAP( dims = 1:10,verbose =F)
  homotypic.prop <- modelHomotypic(tmp@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.05*(ncol(tmp)) ) ## Assuming 5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:10, pN = 0.2, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  tmp=tmp@meta.data
  colnames(tmp)[6]='PANN'
  colnames(tmp)[7]='DF'


  choose.barcode=tmp[tmp$DF=='Singlet',] %>% rownames()
  counts1=counts[,choose.barcode]


saveRDS(counts1,'/data4/huanggy/CTTN/Alldata.removedoublet.rds')




#------cluster CTC step1--------
CTC_counts <- CreateSeuratObject(counts =counts1)
CTC_counts[["percent.mt"]] = PercentageFeatureSet(CTC_counts, pattern = "^[Mm]T-")
CTC_counts[["percent.RPL"]] = PercentageFeatureSet(CTC_counts, pattern = "^RP[SL].+")

CTC_counts <- NormalizeData(object = CTC_counts, normalization.method = "LogNormalize", scale.factor = 10000)
CTC_counts <- FindVariableFeatures(object = CTC_counts, selection.method = "vst", nfeatures = 2000)
CTC_counts <- ScaleData(CTC_counts, vars.to.regress = c("nCount_RNA"))
CTC_counts <- RunPCA(object =CTC_counts, features = VariableFeatures(object = CTC_counts))




message1=rownames(CTC_counts@meta.data) %>% str_split_fixed(';',n=7)
message=rownames(CTC_counts@meta.data) %>% str_split_fixed('_',n=7) 
CTC_counts@meta.data[,'Dataset']=message1[,2]
CTC_counts@meta.data[,'Cancer']=message[,1];CTC_counts@meta.data[,'Patient']=message[,2];CTC_counts@meta.data[,'Tissue']=message[,3]
CTC_counts <- RunHarmony(CTC_counts, group.by.vars = c("Patient"), dims.use = 1:33, verbose = T)
CTC_counts <- FindNeighbors(object = CTC_counts, dims = 1:30, reduction = "harmony")
CTC_counts <- FindClusters(object = CTC_counts, resolution = 0.5)
CTC_counts<- RunUMAP(object = CTC_counts, dims = 1:30, reduction = "harmony")

saveRDS(CTC_counts,'/data4/huanggy/CTTN/Alldata.firstcluster.rds')


#------cluster CTC step2--------
CTC_counts1=subset(CTC_counts1,subset=seurat_clusters%in% c(0,1,3,4))
CTC_counts1 <- NormalizeData(object = CTC_counts1, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
CTC_counts1 <- ScaleData(CTC_counts1, vars.to.regress = c("nCount_RNA"))
CTC_counts1 <- RunPCA(object = CTC_counts1, features = VariableFeatures(object = CTC_counts1))
CTC_counts1 <- RunHarmony(CTC_counts1, group.by.vars = c("Patient"), dims.use = 1:33, verbose = T)
CTC_counts1 <- FindNeighbors(object = CTC_counts1, dims = 1:30, reduction = "harmony")
CTC_counts1 <- FindClusters(object = CTC_counts1, resolution = 0.5)
CTC_counts1<- RunTSNE(CTC_counts1, dims = 1:30, reduction = "harmony")

saveRDS(CTC_counts1,'/data4/huanggy/CTTN/Alldata.Finalcluster.rds')




#----------define CTC subtype-------------
save.pic='/data4/huanggy/CTTN/pic'
markers=wilcoxauc(CTC_counts1,group_by = 'seurat_clusters') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
top2heatmap =markers %>% group_by(group) %>% top_n(n=15,wt=logFC) %>%data.frame() %>% dplyr::arrange(group)
top2heatmap=top2heatmap[!duplicated(top2heatmap$feature),]
marker_heatmap2(top2heatmap,CTC_counts1,save.pic,'seurat_clusters','heatmap.CTC3', levels(as.factor(top2heatmap$group)) ,height=20 )
dev.off()




CTC_counts1@meta.data[CTC_counts1$seurat_clusters%in% c(0), 'recluster' ]='CTC_HBB'
CTC_counts1@meta.data[CTC_counts1$seurat_clusters%in% c(2), 'recluster' ]='CTC_KRT18'
CTC_counts1@meta.data[CTC_counts1$seurat_clusters%in%c(1), 'recluster' ]='CTC_PF4'
CTC_counts1@meta.data[CTC_counts1$seurat_clusters%in%c(3,4), 'recluster' ]='CTC_CTTN'




markers=wilcoxauc(CTC_counts1,group_by = 'recluster') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
top2heatmap =markers %>% group_by(group) %>% top_n(n=15,wt=logFC) %>%data.frame() %>% dplyr::arrange(group)
top2heatmap=top2heatmap[!duplicated(top2heatmap$feature),]
marker_heatmap2(top2heatmap,CTC_counts1,save.pic,'recluster','heatmap.CTC4', levels(as.factor(top2heatmap$group)) ,height=20 )
dev.off()








#---Figure1a left----------------
p=TSNEPlot(CTC_counts1,group.by = 'recluster')+mytheme+scale_color_manual(values = my36colors)
ggsave('/data4/huanggy/CTTN/pic/recluster.pdf',height = 5,width=6,p)

#---Figure1a right--Dot plot--------
Gene.choose <- c('TET2','CTTN','HBB','HBA1','KRT8','KRT18','PPBP','PF4')

p=DotPlot(CTC_counts1,group.by = 'recluster' ,features = Gene.choose ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = '#cacaca', mid = '#ffbf8750', high ='deeppink')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=1, vjust=1,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=60,hjust=1, vjust=1,size = rel(1),color = 'black'))


ggsave('/data4/huanggy/CTTN/pic/Dot.pdf',height = 3,width=4,p)


#---ROE analysis----Figure1b-----------


now.meta=CTC_counts1@meta.data
Roe=Roeforcol(now.meta,'Cancer')
bk = unique(c(seq(0,2, length=100)))
p=pheatmap(Roe,cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = 'black',breaks=bk,fontsize = 12, color=colorRampPalette(c("blue" , "white","orange"),bias=1)(100)) %>% ggplotify::as.ggplot()
ggsave('/data4/huanggy/CTTN/pic/CTC_Roe.pdf',height = 5,width=6,p)


#---FigureS1a---------
p=TSNEPlot(CTC_counts1,group.by = 'Dataset')+mytheme+scale_color_manual(values = my36colors)
ggsave('/data4/huanggy/CTTN/pic/Dataset.pdf',height = 5,width=6,p)

#---FigureS1b---------
p=TSNEPlot(CTC_counts1,group.by = 'Cancer')+mytheme+scale_color_manual(values = my36colors)
ggsave('/data4/huanggy/CTTN/pic/Cancer.pdf',height = 5,width=6,p)




