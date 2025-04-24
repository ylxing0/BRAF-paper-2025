install.packages("tidyverse")
install.packages("readxl")
library(readxl)


### check your data

data1 <- read_excel("Ctrl1.xlsx")
dim(data1) #25239     2
head(data1)

data2 <- read_excel("Ctrl2.xlsx")
dim(data2) #25239     2
head(data2)

identical(rownames(data1),rownames(data2))
#TRUE

data3 <- read_excel("DT1.xlsx")
dim(data3) #25239     2
head(data3)

data4 <- read_excel("DT2.xlsx")
dim(data4) #25239     2
head(data4)

identical(rownames(data1),rownames(data4))
#TRUE

## creat a data matrix
data.mtx<-matrix(data=0,nrow=25239,ncol=4)
dim(data.mtx)
#25239  4
head(data.mtx)

data.mtx[,1]<-as.matrix(data1[,2])
data.mtx[,2]<-as.matrix(data2[,2])
data.mtx[,3]<-as.matrix(data3[,2])
data.mtx[,4]<-as.matrix(data4[,2])

dim(data.mtx)
#25239   4
head(data.mtx)

f1<-c("Ctrl1","Ctrl2","DT1","DT2")
colnames(data.mtx)<-f1
head(data.mtx)

### add more gene information
Geneid.ENSG<-data.frame(data1[,1])
head(Geneid.ENSG)

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("AnnotationDbi")

library(AnnotationDbi)
library(org.Mm.eg.db)

# We can find out which keys and/or columns can be accessed by running
# keytypes(org.Hs.eg.db) or columns(org.Hs.eg.db)

ENTREZID_org <- keys(org.Mm.eg.db, keytype = "ENTREZID")
ENSEMBL_org <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
Symbol_org <- keys(org.Mm.eg.db, keytype = "SYMBOL")

# gene dataframe for org.Mm.eg.db
gene_dataframe_orgDb <- AnnotationDbi::select(org.Mm.eg.db, 
        keys=ENTREZID_org, columns=c("ENSEMBL", "SYMBOL"), keytype="ENTREZID")
colnames(gene_dataframe_orgDb) <- c("Entrez", "Ensembl", "Symbol")

dim(gene_dataframe_orgDb) #73322    3
head(gene_dataframe_orgDb)
data.mtx.inf<-gene_dataframe_orgDb[match(Geneid.ENSG$Geneid,gene_dataframe_orgDb$Symbol),]
dim(data.mtx.inf)  #25239    3
head(data.mtx.inf)
tail(data.mtx.inf)

counts.all<-cbind(data.mtx.inf,data.mtx)
head(counts.all)
dim(counts.all)  #25239    7

## remove the rows with NA
counts.out<-counts.all[complete.cases(counts.all),]
dim(counts.out) # 22978    7
head(counts.out)
write.csv(counts.out,file="Stanford1_dataCountMatrix_mouse.csv")

counts.out <- read.csv("Stanford1_dataCountMatrix_mouse.csv")
head(counts.out[,4:7])

counts.inf<-counts.out[,1:3]
head(counts.inf)
counts<-counts.out[,4:7]
head(counts)
rownames(counts)<-counts.inf$Symbol


#install DESeq2
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("GenomeInfoDb")

library("GenomeInfoDb")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("DESeq2")
#suppressMessages(library(DESeq2))
library("DESeq2")

samples<-data.frame(colnames(counts))
conditions<-c("Ctrl","Ctrl","DT","DT")
samples$conditions<-conditions
samples

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design= ~conditions)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="conditions_DT_vs_Ctrl")


head(res)
summary(res)
#out of 16712 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 856, 5.1%
#LFC < 0 (down)     : 536, 3.2%
#outliers [1]       : 0, 0%
#low counts [2]     : 3176, 19%
#(mean count < 5)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

sum(res$padj < 0.05, na.rm=TRUE)
# 1034


### remove genes without results, add gene information, sort and save 
res.out<-res[complete.cases(res),]
head(res.out)
res.out$minus_log10_adj.pval <- (-log10(res.out$padj))
idx<-match(rownames(res.out),counts.inf$Symbol)
res.out.inf<-cbind(counts.inf[idx,],data.frame(res.out))
dim(res.out.inf) #13536     10

res.sort<-res.out.inf[order(res.out.inf$padj),]
rownames(res.sort)<-NULL
write.csv(res.sort,file="DESeq2_DT_vs_Ctrl_res.csv")


########## MA plot 
pdf("plotMA_mm_results1.pdf", height=5, width=5)
plotMA(res,ylim=c(-6,6),main="DESeq2_DT_vs_Ctrl_res.csv")
dev.off()



#################
pdf("plotCounts_mm_CSC.pdf", height=6, width=15)
par(mfrow=c(2,3))
plotCounts(dds, gene="Nes", intgroup="conditions") #Nes
plotCounts(dds, gene="Sox2", intgroup="conditions") #Sox2
plotCounts(dds, gene="Nanog", intgroup="conditions") #Nanog
plotCounts(dds, gene="Prom1", intgroup="conditions") #Pou5f1
plotCounts(dds, gene="Pou5f1", intgroup="conditions") #Prom1
plotCounts(dds, gene="Cd44", intgroup="conditions") #CD44
dev.off()

###################
pdf("plotCounts_OL.pdf", height=10, width=20)
par(mfrow=c(3,4))
plotCounts(dds, gene="Pdgfra", intgroup="conditions") 
plotCounts(dds, gene="Cspg4", intgroup="conditions") 
plotCounts(dds, gene="Cspg5", intgroup="conditions") 
plotCounts(dds, gene="Olig1", intgroup="conditions") 
plotCounts(dds, gene="Olig2", intgroup="conditions") 
plotCounts(dds, gene="Sox10", intgroup="conditions")
plotCounts(dds, gene="Ncam1", intgroup="conditions") 
plotCounts(dds, gene="Gpr17", intgroup="conditions") 
plotCounts(dds, gene="Galc", intgroup="conditions") 
plotCounts(dds, gene="Ephb1", intgroup="conditions") 
plotCounts(dds, gene="Pcdh7", intgroup="conditions") 
plotCounts(dds, gene="Sox8", intgroup="conditions") 
plotCounts(dds, gene="Slc1a1", intgroup="conditions") 
plotCounts(dds, gene="Plp1", intgroup="conditions") 
plotCounts(dds, gene="Cnp", intgroup="conditions") 
plotCounts(dds, gene="Omg", intgroup="conditions") 
plotCounts(dds, gene="Pcdh9", intgroup="conditions") 
plotCounts(dds, gene="Fabp5", intgroup="conditions") 
plotCounts(dds, gene="Cldn11", intgroup="conditions") 
plotCounts(dds, gene="Mag", intgroup="conditions")
plotCounts(dds, gene="Myrf", intgroup="conditions") 
plotCounts(dds, gene="Mog", intgroup="conditions") 
plotCounts(dds, gene="Mbp", intgroup="conditions")
dev.off()


#################

library(gplots)
library(RColorBrewer)

######### heatmap
install.packages("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
select

head(rownames(res.sort))
sym<-counts.inf$Symbol[match(res.sort$Symbol[1:20],counts.inf$Symbol)]
ntd <- normTransform(dds)
rownames(ntd)<-sym
head(sym)

str(ntd)
vsd <- vst(dds, blind=FALSE)
rownames(vsd)<-sym
head(vsd)

dfm <- as.data.frame(colData(dds)[,c("colnames.counts.","conditions")])
head(dfm)

pheatmap(assay(vsd)[select,][complete.cases(rownames(assay(vsd)[select,])),], 
         cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=dfm)


pheatmap(assay(vsd)[res.sort$Symbol[1:20],], scale="row",
         cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=dfm)

pheatmap(log1p(assay(ntd))[res.sort$Symbol[1:20],], scale="row",
         cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=dfm)

head(counts)
lcpm<-log2((counts+1)/colSums(counts+1)*1000000)

pheatmap(lcpm[res.sort$Symbol[1:20],], scale="row",
         cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=dfm)


library(gplots)
library(RColorBrewer)

pdf("heatmap_mm_top20b.pdf", height=8, width=5)
heatmap.2(as.matrix(lcpm[res.sort$Symbol[1:20],]),main="Top20",
          scale="row", trace="none", key=TRUE,density.info=c("none"),
         distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          dendrogram="row", srtCol=45, cexCol =1, cexRow=0.8,
          margins = c(7, 8), Colv=FALSE,
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

head(counts.out)
dim(counts.out) #22978     7

######## RUN THIS ONE
IFNg_genes <- c("GBP6","FGL2","EIF2AK2","IL2RB","STAT4","IRF1","NMI","UBE2L6","DHX58","USP18",
                "IFITM3","CASP1","PSMB9","NOD1","CIITA","SOD2","B2M","SRI","CXCL10","CXCL11",
                "SECTM1A","IL6","SOCS3","CD74","DDX60","IFIT3","TRIM25","TRIM14","STAT1","H2-AA",
                "XAF1","MX1","OAS2","OAS3","SLC25A28","IFI35","IRF7","H2-AA","RTP4",
                "APOL6","PARP14","C1RB","ZNFX1","SERPING1","RSAD2","HELZ2","CMPK2","STAT2","IFIH1",
                "TNFSF10","HERC6","C1RB","C1S2","H2-Q2","IFI27","NFKBIA",
                "GBP4","PARP12","ISG15","BST2","EPSTI1","PIM1","LGALS3BP","SAMHD1",
                "IFIT2")

cytokines <- c("Ccrl2","Il6ra","Il6","Il1b","Ccl3","Ccl2","Ccl9","Cxcl2","Il1rn","Ccl6","Il17ra","Il22","Cxcl3","Il9", "Il17","Il23","Ifng","Tgfb",
               "Lif","Lifr","Mmp9")
#Above gene list is for mouse genes

library(stringr) 
tmp<-tolower(IFNg_genes)
IFNg_mm_genes <- str_to_title(tmp)
IFNg_mm_genes
temp <- lcpm[cytokines,]
IFNg_mm_genes[!IFNg_mm_genes%in%rownames(temp)]
write.csv(temp,"mouse_IFNg_genes_list.csv")

pdf("heatmap_mm_cytokines.pdf", height=8, width=5)
heatmap.2(as.matrix(temp) ,main="IFNg_genes",
          scale="row", trace="none", key=TRUE,density.info=c("none"),
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          dendrogram="row", srtCol=45, cexCol =1, cexRow=0.5,
          margins = c(7, 8), Colv=FALSE,
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

####################
Tcellactivation <- c("Cd4","Cd8a","Cd3","Foxp3","Tcf7","Cd69","Itgae","Itga1","Cd7","Il2ra",
                     "Tnf","Tbx21","Cd27","Cd28","Ifng","Cxcr6","Cd44","Gzmk","Gzmb","Prf1",
                     "Lck","Klrg1","Pdcd1","Tnfrsf4","Tnfrsf18","Tnfrsf9")
Tcellexhaustion <- c("Pdcd1","Lag3","Ctla4","Havcr2","Tigit","Btla","Eomes","Ccl3","Ccl4","Lgals3","Fas","Socs3")

cytokines <- c("Ccrl2","Il6ra","Il6","Il1b","Ccl3","Ccl2","Ccl9","Cxcl2","Il1rn","Ccl6","Il17ra",
               "Il22","Cxcl3","Il9")



library(stringr)
#Above gene list is for mouse genes



########################
PD1_genes <-c("LCK","CD274","CD3D","PTPN6","CD4","PDCD1","CD247","CD3E",
              "CSK","CD3G","PDCD1LG2","PTPN11")
#Above gene list is for mouse genes
tmp<-tolower(PD1_genes)
PD1_mm_genes <- str_to_title(tmp)
PD1_mm_genes
temp <- lcpm[PDL_mm_genes,]
write.csv(temp,"mouse_PD1_genes_list.csv")

pdf("heatmap_PD1_mm_genes.pdf", height=4, width=3)
heatmap.2(as.matrix(temp),main="PD1_genes",
          scale="row", trace="none", key=TRUE,density.info=c("none"),
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          dendrogram="row", srtCol=45, cexCol =0.8, cexRow=0.8,
          margins = c(7, 8), Colv=FALSE,
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

##################

MHC_genes <-c ("H2-DMa","H2-DMb1","H2-Aa","H2-Ab1")

Interest_genes1 <- c("NES","PROM1","GFAP","CLU","AQP4","CSPG4","PDGFRA","CNP","MBP",
                     "PLP1","PLLP","CLDN11","BTG2","ASPA","OLIG2","EGFR","HES6","ASCL1","DLL1","CD274","LGALS3")



cytokines <- c("Il1b","Ilrn","Il2","Il3","Il4","Il4r","Il5","Il6","Il6r","Il9","Il10","Il12a","Il13","Il17","Il17ra","Il22",
               "Il23a","Il25","Il27","Il29","IL33","Cd274","Il6r")


tmp<-tolower(Interest_genes1)
Interest_mm_genes <- str_to_title(tmp)
Polarity_mm_genes
temp <- lcpm[Candidate,]
temp <- na.omit(temp)
temp <- lcpm[Interest_mm_genes,]


library(circlize)

col_fun = colorRamp2(breaks=c(-2, 0, 2),colors=c("blue","white", "red"))
mat<-t(scale(t(as.matrix(temp))))
mat<-na.omit(mat)


ave.mat<-cbind(( mat[,1]+mat[,2])/2,(mat[,3]+mat[,4])/2)

pdf("heatmap_candidate_glial_genes_ave.pdf", height=5, width=3)
Heatmap(ave.mat, name = "Expression",  
        cluster_columns = TRUE,
        show_column_dend = TRUE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 3),
        column_gap = unit(0.5, "mm"),
        cluster_rows = TRUE,
        show_row_dend = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 6),
        column_title_rot = 45,
        show_column_names = TRUE,
        use_raster = FALSE,
        raster_quality = 20,
        width = unit(4, "cm"), height = unit(7, "cm"))
dev.off()



pdf("Stripchat_Interest_mm_genes1.pdf", height=6,width=4)
nice.col <- brewer.pal(6,name="Set1")
par(mfrow=c(3,3),mar=c(5,5,2,1))

ng2display<-6
#g.list1<-topTab1$ID[1:ng2display]
for (i in 1:ng2display) { 
  stripchart(as.matrix(lcpm)[rownames(lcpm)%in%Interest_mm_genes[i],]~samples$conditions,vertical=TRUE,las=2,cex.axis=1,pch=16,cex=1.3,
             col=mycol,method="jitter",ylab="expression",main=rownames(temp)[i])
}

dev.off()
#########


library(stringr)
library(gplots)
library(RColorBrewer)
##################
#######

pdf("Galectin3_expression_in_2341.pdf", height=3, width=2)
stripchart(as.matrix(lcpm)[rownames(lcpm)%in%c("Lgals3"),]~samples$conditions,vertical=TRUE,las=2,cex.axis=1,pch=16,cex=1.3,
           col=mycol,method="jitter",ylab="expression",main=c("Galectin-3"))
dev.off()

###################

####################################
#### using edgeR com function # edgeR is problematic, unable to install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("impute")
BiocManager::install("edgeR",force = TRUE)


#BiocManager::install("genefilter")

dyn.load('/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/edgeR/libs/edgeR.so')
dyn.load('/opt/R/arm64/lib/libgfortran.5.dylib')
dyn.load('/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/edgeR')
dyn.load('edgeR.so')

#library(impute)
library("edgeR")

lcmp4heatmap<-cpm(counts.out[,4:7],normalized.lib.sizes = TRUE,log=TRUE)
rownames(lcmp4heatmap)<-counts.out[,4]
colnames(lcmp4heatmap)<-gsub("SU-aGBM5_","",samples$colnames.counts.)

pheatmap(lcmp4heatmap[res.sort$Symbol[1:20],],main="top20 DT24_vs_DMSO",
         clustering_method="ward.D2",angle_col=45)
pheatmap(lcmp4heatmap[res2.sort$HGNC[1:20],],main="top20 DT48_vs_DMSO",
         clustering_method="ward.D2",angle_col=45)


##################
#################
#################
library(clusterProfiler) 
library(org.Mm.eg.db) 
library(org.Hs.eg.db) 
library(Seurat) 
library(dplyr)
library(tibble) 
library(GOSemSim) 
library(enrichplot) 
library(AnnotationDbi)


####GO analysis 
mm_genes <- read.csv("DESeq2_DT_vs_Ctrl_res.csv", header=TRUE)
head(mm_genes)
mm_genes_up <-subset(mm_genes,subset=mm_genes$log2FoldChange>0)
mm_genes_down <-subset(mm_genes,subset=mm_genes$log2FoldChange<0)

mm_genes_up1<-mm_genes_up$Entrez
mm_genes_down1<-mm_genes_down$Entrez

mm_genes_up1<-list(X0=mm_genes_up1)
mm_genes_down1<-list(mm_genes_down1)

lfc_up <-mm_genes_up$log2FoldChange
names(lfc_up) <-mm_genes_up$Entrez
lcf_up<-na.omit(lfc_up)
lcf_up<-sort(lcf_up,decreasing=TRUE)
head(lcf_up)

lfc_down <-mm_genes_down$log2FoldChange
names(lfc_down) <-mm_genes_down$Entrez
lcf_down<-na.omit(lfc_down)
lcf_down<-sort(lcf_down,decreasing=TRUE)
head(lfc_down)
tail(lfc_down)

lcf <-c(lcf_up,lcf_down)
head(lcf)
tail(lcf)
lcf<-na.omit(lcf)

tail(mm_genes_down1)
Combined_genelist<-cbind(mm_genes_up1,mm_genes_down1)
names(Combined_genelist)<-c("X0","X1")

x<-as.matrix(Combined_genelist)
Combined_genelist<-tapply(x,rep(1:ncol(x),each=nrow(x)),function(i)i)
str(Combined_genelist)
class(Combined_genelist)
head(Combined_genelist)
Combined_genelist <-na.omit(Combined_genelist)

lapply(Combined_genelist, head)
ont<-c("MF","CC","BP")
j=1
ck_up<-compareCluster(geneClusters=Combined_genelist,fun="enrichGO", 
                      pvalueCutoff=0.05, OrgDb=org.Mm.eg.db, ont=ont[j])
head(as.data.frame(ck_up))
pdf("Dotplot_mm_genes_GO_MF.pdf", width=5, height=6)
dotplot(ck_up, showCategory=10, font.size=8)
dev.off()


####
####

ont<-c("MF","CC","BP")
j=3
ego1 <- enrichGO(gene     = Combined_genelist$`1`,
                 OrgDb    = org.Mm.eg.db,
                 ont      = ont[j],
                 pAdjustMethod="BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = TRUE)


head(ego1)
#pdf(paste0("Heatplot_mmBP_up.pdf"),height=5,width=15)
#heatplot(ego1)
#dev.off()

pdf(paste0("Cnetplot_mm_BP_up.pdf"),height=20,width=25)
#cnetplot(ego1,categorySize="p.adjust")
cnetplot(ego1,categorySize="p.adjust", colorEdge=TRUE, 
         foldChange=lcf_up,showCategory=5)
dev.off()
#######


d<-godata('org.Mm.eg.db', ont="BP")
ego2<-pairwise_termsim(ego1,method="Wang",semData=d)
ck1_up<-pairwise_termsim(ck_up,method="Wang",semData=d)

pdf("emapplot_mm_BPa_up.pdf",height=20, width=20)
emapplot(ego2)
dev.off()

pdf("emapplot_mm_BPb_up.pdf",height=10, width=10)
emapplot(ck1_up)
dev.off()

##############
ont<-c("MF","CC","BP")
j=3
ego3 <- enrichGO(gene     = Combined_genelist$`2`,
                 OrgDb    = org.Mm.eg.db,
                 ont      = ont[j],
                 pAdjustMethod="BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = TRUE)
#######
pdf(paste0("Cnetplot_mm_BP_down1.pdf"),height=10,width=15)
#cnetplot(ego1,categorySize="p.adjust")
cnetplot(ego3,categorySize="p.adjust", colorEdge=TRUE, foldChange=lcf_down,showCategory=5)
dev.off()
#######


d<-godata('org.Mm.eg.db', ont="BP")
ego4<-pairwise_termsim(ego3,method="Wang",semData=d)
ck2_up<-pairwise_termsim(ck_up,method="Wang",semData=d)

pdf("emapplot_mm_BPa_down.pdf",height=12, width=12)
emapplot(ego4)
dev.off()

pdf("emapplot_mm_BPb_down.pdf",height=8, width=8)
emapplot(ck2_up)
dev.off()



#######
kk <- enrichKEGG(gene         = Combined_genelist$`1`,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
head(kk)

pdf("Dotplot_mm_genes_KK_up.pdf", width=7, height=7)
dotplot(kk)
dev.off()

########
kk_down <- enrichKEGG(gene         = Combined_genelist$`2`,
                      organism     = 'mmu',
                      pvalueCutoff = 0.05)
head(kk)

pdf("Dotplot_mm_genes_KK_down.pdf", width=6, height=5)
dotplot(kk_down)
dev.off()

kk<- setReadable(kk_down,OrgDb=org.Mm.eg.db,keyType="ENTREZID")

pdf("Cnetplot_ms_genes_Endocytosis", width=20, width=20)
cnetplot(kk, showCategory=c("Endocytosis"), colorEdge=TRUE, foldChange=lcf_down)
dev.off()



####### REACTOME
#install.packages("reactome.db")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("reactome.db", force=TRUE)

#f (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("ReactomePA")

###### REACTOME up
library(ReactomePA)
x <- enrichPathway(gene=Combined_genelist$`1`,pvalueCutoff=0.05, readable=T, organism="mouse")
head(as.data.frame(x))

pdf("Barplot_mm_genes_Reactome_up.pdf", width=4, height=4)
barplot(x, showCategory=10, font.size =6)
dev.off()

pdf("Dotplot_mm_genes_Reactome_up.pdf", width=7, height=8)
dotplot(x, showCategory=10)
dev.off()


pdf("Cnetplot_mm_genes_Reactome_up.pdf", width=10, height=10)
cnetplot(x, categorySize="pvalue", foldChange=lcf_up)
dev.off()

pdf("Treeplot_mm_genes_Reactome_up.pdf", width=20, height=10)
treeplot(pairwise_termsim(x))
dev.off()

res <- compareCluster(Combined_genelist, fun="enrichPathway", organism="mouse")
pdf("Dotplot_mm_genes_Reactome_Up_down.pdf", width=5, height=6)
dotplot(res)
dev.off()

########### REACTOME Down
x1 <- enrichPathway(gene=Combined_genelist$`2`,pvalueCutoff=0.05, readable=T, organism="mouse")
head(as.data.frame(x1))

pdf("Barplot_mm_genes_Reactome_down.pdf", width=4, height=4)
barplot(x1, showCategory=10, font.size =6)
dev.off()

pdf("Dotplot_mm_genes_Reactome_down.pdf", width=7, height=8)
dotplot(x1, showCategory=10)
dev.off()


pdf("Cnetplot_mm_genes_Reactome_down.pdf", width=15, height=15)
cnetplot(x1, categorySize="pvalue", foldChange=lcf_down)
dev.off()

pdf("Treeplot_mm_genes_Reactome_down.pdf", width=20, height=10)
treeplot(pairwise_termsim(x1))
dev.off()


############
y <- gsePathway(lcf, nPerm=10000,
                pvalueCutoff=0.2,
                pAdjustMethod="BH", verbose=FALSE, organism = "mouse")
yy<- setReadable(y,OrgDb=org.Mm.eg.db,keyType="ENTREZID")
res <- as.data.frame(y)
head(y)

y2<-pairwise_termsim(y,method="JC")

pdf("emapplot_mm_GSEA_REACTOME.pdf",height=10, width=20)
emapplot(y2, color="pvalue")
dev.off()

pdf("Dotplot_mm_GSEA_REACTOME_paper.pdf",height=3, width=3)
dotplot(yy, showCategory=10, font.size=5)
dev.off()

write.csv(yy,file="mm_GSEA_results_REACTOME.csv")


pdf("Gseaplot_PD-1_signaling.pdf", height=3, width=4)
gseaplot(yy, geneSetID = "R-MMU-389948")
gseaplot2(yy, geneSetID= "R-MMU-389948")
dev.off()

pdf("Gseaplot_top10a.pdf", height=18, width=15)
gseaplot2(yy, geneSetID= 1:20, pvalue_table=TRUE, base_size = 20)
dev.off()

?gseaplot2


pdf("Viewpathway_mm_genes_GSEA_PD-1_signaling1.pdf", width=5, height=4)
viewPathway("PD-1 signaling", readable=TRUE, foldChange=lcf[unique(names(lcf))], organism="mouse") +
  scale_color_continuous(low="blue", high="red", name = "fold change", na.value = "#E5C494")
dev.off()

pdf("Viewpathway_mm_genes_GSEA_ER-Phagosome pathway.pdf", width=5, height=4)
viewPathway("ER-Phagosome pathway", readable=TRUE, foldChange=lcf[unique(names(lcf))], organism="mouse")
dev.off()

pdf("Ridgeplot_mm_GSEA_reactome.pdf", width=6, height=20)
ridgeplot(yy)
dev.off()

pdf("Treeplot_mm_GSEA_reactome.pdf", width=20, height=10)
treeplot(pairwise_termsim(yy))
dev.off()


?ridgeplot
?cnetplot
?viewPathway

############## Don't run
organism="org.Hs.eg.db"
gse <- gseGO(geneList=lcf[abs(lcf)>1.5], 
             ont ="ALL", 
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

dotplot(gse, showCategory=10)
emapplot(pairwise_termsim(gse), showCategory = 10)

###################
kegg_organism = "mmu"
kk2 <- gseKEGG(geneList     = lcf[abs(lcf)>1.5], 
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways")


############

install.packages("msigdbr")
library(msigdbr)
msigdbr_show_species()
H <- msigdbr(species="Mus musculus", category="H") %>%
  dplyr::select(gs_name,entrez_gene)

en1<-enricher(names(lcf_up), TERM2GENE = H)
en1<-setReadable(en1,OrgDb=org.Mm.eg.db,keyType = "ENTREZID")

pdf("mm_hallmark_results_a.pdf", height=3, width=6)
barplot(en1,showCategory = 20, font.size =5)
dotplot(en1,showCategory = 20, font.size=5,size=NULL)
dev.off()

?dotplot

pdf("mm_hallmark_results_b.pdf", height=5, width=16)
treeplot(pairwise_termsim(en1), font.size=4)
dev.off()

####################




