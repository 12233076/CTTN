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
save.data='/data4/huanggy/CTTN/Inhouse_CTC_singleron_validation/data/'
save.pic='/data4/huanggy/CTTN/Inhouse_CTC_singleron_validation/pic/'


#-----loading data------

MEL167.PM <- readRDS("/data5/huanggy/CTC_database_SC/MEL167.PM.rds")

#----cluster CTCs celline-----step1-------
MEL167.Cellline=subset(MEL167.PM,subset=Tissue%in% c('Cellline'))
MEL167.Cellline <- NormalizeData(object = MEL167.Cellline, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
MEL167.Cellline <- ScaleData(MEL167.Cellline, vars.to.regress = c("nCount_RNA"))
MEL167.Cellline <- RunPCA(object = MEL167.Cellline, features = VariableFeatures(object = MEL167.Cellline))
MEL167.Cellline <- FindNeighbors(MEL167.Cellline,dims = 1:35, reduction = "harmony") %>% FindClusters(resolution = 0.5)
MEL167.Cellline <- RunUMAP(object = MEL167.Cellline, dims = 1:30, reduction = "harmony")


#----cluster CTCs celline-----step2-------
MEL167.Cellline_CTTN=subset(MEL167.Cellline,subset=Name%in% c('high'))
MEL167.Cellline_CTTN <- NormalizeData(object = MEL167.Cellline_CTTN, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
MEL167.Cellline_CTTN <- ScaleData(MEL167.Cellline_CTTN, vars.to.regress = c("nCount_RNA"))
MEL167.Cellline_CTTN <- RunPCA(object = MEL167.Cellline_CTTN, features = VariableFeatures(object = MEL167.Cellline))
MEL167.Cellline_CTTN <- FindNeighbors(MEL167.Cellline_CTTN,dims = 1:35, reduction = "harmony") %>% FindClusters(resolution = 0.5)
MEL167.Cellline_CTTN <- RunUMAP(object = MEL167.Cellline_CTTN, dims = 1:30, reduction = "harmony")


#Base on the epxression of CTTN of secondary CTCs and seperate into 2 groups(Top 25% and Bottom 25%)------------------
nouse1 <- MEL167.Cellline_CTTN@assays$RNA@data%>%t() %>%data.frame()
MEL167.Cellline_CTTN@meta.data$CTTN <- nouse1[rownames(MEL167.Cellline_CTTN@meta.data),'CTTN']
MEL167.Cellline_CTTN@meta.data$Name <- 'Neu'
quantile(MEL167.Cellline_CTTN@meta.data$CTTN,0.25)
MEL167.Cellline_CTTN@meta.data[MEL167.Cellline_CTTN@meta.data$CTTN>= 2.393589,'Name']='high'
MEL167.Cellline_CTTN@meta.data[MEL167.Cellline_CTTN@meta.data$CTTN < 0.6629637,'Name']='low'



#----cluster CTCs celline-----step3-------
MEL167.Cellline_CTTN1 <-subset(MEL167.Cellline_CTTN,subset=Name%in% c('high','low'))
MEL167.Cellline_CTTN1 <- NormalizeData(object = MEL167.Cellline_CTTN1, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
MEL167.Cellline_CTTN1 <- ScaleData(MEL167.Cellline_CTTN1, vars.to.regress = c("nCount_RNA"))
MEL167.Cellline_CTTN1 <- RunPCA(object = MEL167.Cellline_CTTN1, features = VariableFeatures(object = MEL167.Cellline))
MEL167.Cellline_CTTN1 <- FindNeighbors(MEL167.Cellline_CTTN1,dims = 1:35, reduction = "harmony") %>% FindClusters(resolution = 0.5)
MEL167.Cellline_CTTN1 <- RunUMAP(object = MEL167.Cellline_CTTN1, dims = 1:30, reduction = "harmony")





CTC_counts1=MEL167.Cellline_CTTN1


# cell cycle score---------------
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

CTC_counts1 <- CellCycleScoring(CTC_counts1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#Define potential Dormancy CTCs-----------

# Plot the two scores, and draw two lines to indicate the thresholds Score_S + Score_G2M > -0.2
# (let's consider these cycling) and Score_S + Score_G2M < -0.5 (let's consider these G0)

CTC_counts1$CycleSum <- CTC_counts1@meta.data$S.Score+CTC_counts1@meta.data$G2M.Score
CTC_counts1$Cycle <- 'G1/G0'
CTC_counts1@meta.data[CTC_counts1@meta.data$CycleSum>= -.2,'Cycle']='cycling'
CTC_counts1@meta.data[CTC_counts1@meta.data$CycleSum < -.5,'Cycle']='G0'

#--------------------Figure S2a-----------
p=ggplot( CTC_counts1@meta.data ) +
  geom_point( aes( x=S.Score, y=G2M.Score, col=Name ) ) +xlim(-1,1)+ylim(-1,1)+
  coord_fixed() +
  geom_abline( intercept = c(-1,-.5), slope=-1, col="gray" )+scale_color_manual(values=my36colors)
ggsave('/data4/huanggy/CTTN/Inhouse_CTC_singleron_validation/pic/cellcycle-tissue.pdf',height = 4,width=4.5,p)


#--------------------Figure S2b-----------
tmp.pic=CTC_counts1@meta.data %>% data.table() %$% .[,.N,.(Cycle,Name)] %>% data.frame() %>% ddply('Cycle',transform,per=N*100/sum(N))
p=ggplot(tmp.pic,aes(x=Cycle,y=per,fill=Name))+geom_bar(stat='identity')+mytheme+scale_fill_manual(values = c(my36colors))+
  theme(axis.text.x = element_text(angle=30,hjust=1, vjust=1,size = rel(1.2)))


ggsave('/data4/huanggy/CTTN/Inhouse_CTC_singleron_validation/pic/cellcycle_bar-CTTN.pdf',height = 4,width=4,p)


#--------------------Figure S2c-----------
StackedVlnPlot(CTC_counts1,c('CDKN1A','MMP14','ICAM1','IGFBP3','MIF4GD','VEGFA','CD55','PLAT','SERPINE2','GDF15'), 'Name',color=my36colors,pt.size=1)






#irgesa------
markers=list()

kegmt<-read.gmt('/data4/huanggy/GSEA_gmt/c5.go.v7.4.symbols.gmt')
kegmt[,1]=kegmt[,1] %>% str_remove('GO_')  %>% tolower()%>% capitalize()

#Cellular_senescence
markers$Cellular_senescence=kegmt[kegmt$term=='Gobp_cellular_senescence',2] 
#cellcycle_arrest
markers$cellcycle_arrest=kegmt[kegmt$term=='Gobp_cell_cycle_arrest',2] 

kegmt<-read.gmt('/data4/huanggy/GSEA_gmt/c2.cp.kegg.v7.4.symbols.gmt')
kegmt[,1]=kegmt[,1] %>% str_remove('KEGG_')  %>% tolower()%>% capitalize()

#ROS
markers$Reactive_oxygen_species_pathway= kegmt[kegmt$term=='Reactive_oxygen_species_pathway',2] 
#P53
markers$P53_pathway= kegmt[kegmt$term=='P53_pathway',2] 
#Oxidative_phosphorylation
markers$Oxidative_phosphorylation= kegmt[kegmt$term=='Oxidative_phosphorylation',2]
#senmayo
senmayo<-read.table('/data4/huanggy/jinayang_SC/MEL167/senmayo_geneset.csv',header = T,sep = ',',row.names = 1)
senmayo$name<-rownames(senmayo)
#CTTNKD-induced SASP
senmayo_upregulate <- c('IL1B','PLAUR','VGF','C3','TNFRSF1A','GDF15','MIF4GD','MMP3','CCL2','FAS','CXCL8','CD55','DKK1','IGFBP3','PLAT','SERPINE2','PGF','ICAM1','MMP14','SERPINE1','CXCL1','VEGFA','CD96','ANGPTL4','CTSB')

markers$senmayo<-senmayo$name
markers$senmayo_upregulate<-senmayo_upregulate
#Mitochondrial.ETC
ETC <- read.gmt('/data4/huanggy/GSEA_gmt/WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA.v2023.1.Hs.gmt')
markers$mitochondrial_ETC = ETC$gene



count.one.type.pathway.CTC <- irGSEA.score(object =CTC_counts1, assay = "RNA", slot = "data", seeds = 123, ncores = 10,
                                            min.cells = 3, min.feature = 0,custom = T, geneset = markers, msigdb = F,
                                            species = "Homo sapiens", category = "H",  subcategory = NULL, geneid = "symbol",
                                            method = c("AUCell","ssgsea"),aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                            kcdf = 'Gaussian')

#Boxplot----------FigureS 2d--------------------
nouse <-count.one.type.pathway.CTC@assays$AUCell@data%>%t() %>% data.frame()

genes <- colnames(nouse)


CTC_counts1@meta.data[,genes] <- nouse[rownames(CTC_counts1@meta.data),genes]


df = data.frame(t(count.one.type.pathway.CTC@assays$AUCell$scale.data))
df = df[names(CTC_counts1$Cycle),]

df$Name <- CTC_counts1$Cycle


df_melt <- melt(df,ID='Name')
p=ggplot(df_melt,aes(x=variable,y=value,fill=Name))+geom_boxplot()+mytheme+ylim(-0.1,0.7)+scale_fill_manual(values=my36colors)+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Name),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=1, vjust=1,size = rel(1)))
ggsave('/data5/huanggy/CTTN_project/SC_Validation/SC_G0.pdf',p,height = 4,width=4)



