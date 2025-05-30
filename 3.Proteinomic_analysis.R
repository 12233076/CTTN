library(stringr)
library(DESeq2)
library(ggrepel)
library(ggplot2)
library(edgeR)
library(limma)
library(dplyr)
library(clusterProfiler)
has_kegg<-download_KEGG('hsa')
library(DOSE)
library(topGO)
library(pathview) 
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichplot)
library(RColorBrewer)
library(GseaVis)
library(VennDiagram)
library(gplots)
library(Hmisc)
library(KEGGREST)
library(plyr)
library(parallel)
library(corrplot)
library(pheatmap)
library(ggvenn)

set.seed(1234)

my36colors <-c('#D6E7A3',"orange1", 'lightblue','#7142AC',"darkcyan","royalblue1","red3",'#53A85F',"deeppink",
               "mediumvioletred","gold","darkorange2", "tan2","darkorange2","darkorchid","chocolate4","darkred","lightskyblue","gold1")

mytheme <- theme_bw() + 
  theme(plot.title=element_text(size=rel(2),hjust=0.5),
        axis.title=element_text(size=rel(1)),
        axis.text.x = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2),colour = 'black'),
        axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2),colour = 'black'),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        panel.border=element_rect(color="black",size=1),
        axis.line=element_line(color="black",size=0.5))

main.path='/data5/huanggy/jianyang_CTTN/'
save.pic='/data5/huanggy/jianyang_CTTN/pic/'
save.data='/data5/huanggy/jianyang_CTTN/'

#----matrix loading-------

protein <- read.csv('/data5/huanggy/jianyang_CTTN/Protein/protein_expression.csv',header = T,sep = ',')
protein[protein== ""] <- NA
protein <- na.omit(protein)
protein <- protein[!duplicated(protein$Genes),]

rownames(protein) <-protein$Genes
protein$Genes =NULL


protein=protein[,c(7:9,4:6)]
gene_select=data.frame(row.names = colnames(protein),FDPS=protein['FDPS',] %>% as.numeric())
gene_select$Group=c('CTRL','CTRL','CTRL','CTTNko','CTTNko','CTTNko')
gene_select$name=rownames(gene_select)

#----DEG analysis---------
y <- DGEList(counts=as.matrix(protein),group=gene_select$Group )
y <- y[rowSums(edgeR::cpm(y)>1) >= 2, , keep.lib.sizes=FALSE] %>% calcNormFactors() %>% estimateDisp( model.matrix(~gene_select$Group))
et <- exactTest(y)
res <- et$table %>% dplyr::mutate( name=rownames(.) ) %>% dplyr::mutate(logpvalue= -log10(.$PValue) ) %>%  { .[! .$name %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]} %>%
  { .[is.infinite(.$logpvalue),'logpvalue'] = max(.$logpvalue[.$logpvalue!=max(.$logpvalue)]) ;. } %>%  {.[.$logpvalue>300,'logpvalue']=300 ;.}
res = res %>% {.[.$logFC>0.58&.$ logpvalue>1.3,'color']='CTTNkoupregulate';.}%>% {.[.$logFC< -0.58&.$logpvalue>1.3,'color']='CTTNkodownregulate';. } %>% {.[is.na(.$color) ,'color']= 'No Sig';.}

#----protein DEGs for downstream analysis----------
gene_diff_CTRL_CTTNko_protein<-res[,c('logFC','logpvalue','color','name')]


#-heatmap--Figure3a------------
kegmt<-read.gmt('/data5/huanggy/TCGA_database/data/GSEA_gmt/h.all.v7.0.symbols.gmt')
kegmt[,1]=capitalize(tolower(str_remove(kegmt[,1],'HALLMARK_')))
markers=list()
#P53
markers$P53_pathway= kegmt[kegmt$term=='P53_pathway',2] 

heatmap_protein<-protein[markers$P53_pathway,]
heatmap_protein<-na.omit(heatmap_protein)
heatmap_protein<-cbind(heatmap_protein,Total=rowSums(heatmap_protein))
heatmap_protein<-heatmap_protein[order(heatmap_protein$Total),]
heatmap_protein<-heatmap_protein[!is.na(heatmap_protein$Total),]
heatmap_protein[which(heatmap_protein$Total > 0),'Group']<-'Right'
heatmap_protein[which(heatmap_protein$Total <= 0),'Group']<-'Left'
heatmap_protein<-subset(heatmap_protein,subset= Group%in%c('Right'))
heatmap_protein$Total<-NULL
heatmap_protein$Group<-NULL

bk = unique(c(seq(-1,1, length=100)))
pdf(paste0(save.data,'bulk_protein.heatmap_P53_pathway.pdf'),height = 30,width = 6)
pheatmap(heatmap_protein,scale='row',breaks = bk,cluster_cols = F,show_colnames =T,show_rownames = T,color=colorRampPalette(c("blue", "black", "yellow"))(100))
dev.off()

##----enrichment analysis ---Figure3b------
tmp<-gene_diff_CTRL_CTTNko_protein
kegmt<-read.gmt('/data5/huanggy/TCGA_database/data/GSEA_gmt/h.all.v7.0.symbols.gmt')
kegmt[,1]=capitalize(tolower(str_remove(kegmt[,1],'HALLMARK_')))

hark<-GSEA(tmp$logFC %>% { names(.)=tmp$name;.} %>% sort(decreasing = T),TERM2GENE = kegmt) 

##----enrichment analysis  hallmark50---Figure S4b----
p=ggplot(hark@result,aes(x=NES,y=reorder(Description,NES)))+geom_segment(aes(yend=Description),xend=0,colour='grey50')+mytheme+geom_point(size=5,aes(colour = p.adjust))+scale_color_gradientn(colors = c('red','#ffbf8750','#cacaca'))+theme(text=element_text(family="sans",face = 'bold'))+theme(axis.text.y = element_text(hjust = 1))

ggsave('/data5/huanggy/CTTN_project/jianyang_CTTN/protein_Hallmark50.pdf',height = 5,width=7,p)


p=gseaplot2(hark,geneSetID = 'P53_pathway',pvalue_table=TRUE)

ggsave(paste0(save.pic,'Hallmark_P53_pathway.pdf'),width = 6,height =6,p)


#melanoma 500 patient public protein database ---Figure3d -------

j='Q14247'#CTTN
j='P04637'#TP53


protein <- read.csv('/data5/huanggy/CTTN_project/jianyang_CTTN/Protein/public_protein_database.csv',header = T,sep = ',')

protein$Peptide <- str_replace_all(protein$Peptide,'-','_')
rownames(protein) <-protein$Peptide
annotation <- protein[,c(1:3)]%>%data.frame()
protein$X =NULL
protein$Peptide =NULL



protein1<- protein
protein1$Name=NULL
protein1<- protein1%>%t()%>%data.frame()



protein.COR <- protein1[!is.na(protein1$Q14247),]
protein.COR <- protein.COR[!is.na(protein.COR$P04637),]


p=ggplot(protein.COR,aes(x=Q14247,y=P04637))+geom_point(size=3,shape=19,color='lightblue',alpha=0.5)+mytheme+geom_smooth(method ='lm',color='mediumvioletred')+stat_cor(method='spearman')+scale_color_manual(values=my36colors)





#melanoma 500 patient public protein database ---Figure S4h -------
'Q14247'#CTTN
'P38936'#CDKN1A

protein1<- protein
protein1$Name=NULL
protein1<- protein1%>%t()%>%data.frame()

protein.COR <- protein1[!is.na(protein1$Q14247),]
protein.COR <- protein.COR[!is.na(protein.COR$P38936),]


p=ggplot(protein.COR,aes(x=Q14247,y=P38936))+geom_point(size=3,shape=19,color='lightblue',alpha=0.5)+mytheme+geom_smooth(method ='lm',color='mediumvioletred')+stat_cor(method='spearman')+scale_color_manual(values=my36colors)




























