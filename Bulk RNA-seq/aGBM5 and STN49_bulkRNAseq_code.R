install.packages("tidyverse")
install.packages("readxl")

library(readxl)

### learn to read excel file
readxl_example()
xls_example <- readxl_example("datasets.xls")
read_excel(xls_example)
excel_sheets(xls_example)
read_excel(xls_example, sheet = "chickwts")


### pre-processing of aGBM5 files

dir("tmp/Manuscript/")
f.names<-dir("tmp/Manuscript/")
f.names
#[1] "SU-aGBM5_DMSO_1.xlsx" "SU-aGBM5_DMSO_2.xlsx" "SU-aGBM5_DT24_1.xlsx"
#[4] "SU-aGBM5_DT24_2.xlsx" "SU-aGBM5_DT48_1.xlsx" "SU-aGBM5_DT48_2.xlsx"

f.names[1]
#"SU-aGBM5_DMSO_1.xlsx"

data1<-read_excel(paste0("tmp/Manuscript/",f.names[1]))
dim(data1) #57773     2
head(data1)

data1<-read_excel(paste0("tmp/Manuscript/",f.names[1]),range = cell_rows(1:57773))
dim(data1) #57772     2
head(data1)

data2<-read_excel(paste0("tmp/Manuscript/",f.names[2]),range = cell_rows(1:57773))
dim(data2) #57772     2

identical(rownames(data1),rownames(data2))
#TRUE

data6<-read_excel(paste0("tmp/Manuscript/",f.names[6]),range = cell_rows(1:57773))
dim(data6) #57772     2

identical(rownames(data1),rownames(data6))
#TRUE

## creat a data matrix
data.mtx<-matrix(data=0,nrow=57772,ncol=6)

dim(data.mtx)
#57772     6
head(data.mtx)

data.mtx[,1]<-as.matrix(data1[,2])
data.mtx[,1]
for (i in 2:6) {
  data<-read_excel(paste0("tmp/Manuscript/",f.names[i]),range = cell_rows(1:57773))
  data.mtx[,i]<-as.matrix(data[,2])
}

dim(data.mtx)
#57772     6
head(data.mtx)

f1<-gsub("SU-aGBM5_","",f.names)
f2<-gsub(".xlsx","",f.names)
f2
colnames(data.mtx)<-f2
head(data.mtx)

### add more gene information
Geneid.ENSG<-data.frame(data1[,1])

library(AnnotationDbi)
library(org.Hs.eg.db)

# We can find out which keys and/or columns can be accessed by running
# keytypes(org.Hs.eg.db) or columns(org.Hs.eg.db)

ENTREZID_org <- keys(org.Hs.eg.db, keytype = "ENTREZID")
ENSEMBL_org <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
HGNC_org <- keys(org.Hs.eg.db, keytype = "SYMBOL")

# gene dataframe for org.Hs.eg.db
gene_dataframe_orgDb <- AnnotationDbi::select(org.Hs.eg.db, keys=ENTREZID_org, columns=c("ENSEMBL", "SYMBOL"), keytype="ENTREZID")
colnames(gene_dataframe_orgDb) <- c("Entrez", "Ensembl", "HGNC")

dim(gene_dataframe_orgDb) #69580     3
head(gene_dataframe_orgDb)
data.mtx.inf<-gene_dataframe_orgDb[match(Geneid.ENSG$Geneid,gene_dataframe_orgDb$Ensembl),]
dim(data.mtx.inf)  #57772     3

counts.all<-cbind(data.mtx.inf,data.mtx)
head(counts.all)
dim(counts.all)  #57772     9

## remove the rows with NA
counts.out<-counts.all[complete.cases(counts.all),]
dim(counts.out) # 33158     9
head(counts.out)
write.csv(counts.out,file="Stanford1_dataCountMatrix.csv")

########### RUN THIS CODE
counts.out <- read.csv("Stanford1_dataCountMatrix.csv")
counts.out <- counts.out[!duplicated(counts.out$HGNC),]
head(counts.out[,6:11])

counts.inf<-counts.out[,1:5]
counts<-counts.out[,6:11]
rownames(counts)<-counts.inf$HGNC
head(counts)

class(counts)
counts[,1:6]<-lapply(counts,as.numeric)
counts[is.na(counts)]<-0

lcpm<-log2((counts+1)/colSums(counts+1)*1000000)
head(lcpm)

######################
#install DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDb")

library("GenomeInfoDb")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
suppressMessages(library(DESeq2))


samples<-data.frame(colnames(counts))
conditions<-gsub("SU-aGBM5_","",colnames(counts))
samples$conditions<-substr(conditions,start = 10,stop =13)
samples

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design= ~conditions)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res1 <- results(dds, name="conditions_aGBM5_DT24_vs_aGBM5_DMSO")
res2 <- results(dds, name="conditions_aGBM5_DT48_vs_aGBM5_DMSO")

res<-results(dds, contrast=c("conditions","aGBM5_DT24","aGBM5_DMSO"))

hist(res1$padj, 
     col="grey", border="white", xlab="", ylab="", main="frequencies of adj. p-values\n(all genes)")

head(res1)
summary(res1)
#out of 20003 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2587, 13%
#LFC < 0 (down)     : 2832, 14%
#outliers [1]       : 0, 0%
#low counts [2]     : 5990, 30%
#(mean count < 6)
sum(res1$padj < 0.05, na.rm=TRUE)
# 4682


head(res2)
summary(res2.sort)
#out of 20003 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2990, 15%
#LFC < 0 (down)     : 3342, 17%
#outliers [1]       : 0, 0%
#low counts [2]     : 5241, 26%
#(mean count < 3)
sum(res2$padj < 0.05, na.rm=TRUE)
# 5576

### remove genes without results, add gene information, sort and save 
res1.out<-res1[complete.cases(res1),]
head(res1.out)
res1.out$minus_log10_adj.pval <- (-log10(res1.out$padj))
idx<-match(rownames(res1.out),counts.inf$Ensembl)
head(res1.out.inf)
res1.out.inf<-cbind(counts.inf[idx,],data.frame(res1.out))
dim(res1.out.inf) #14013     11

res1.sort<-res1.out.inf[order(res1.out.inf$padj),]
rownames(res1.sort)<-NULL
write.csv(res1.sort,file="DESeq2_DT24_vs_DMSO_res1.csv")

res1.sort <- read.csv("DESeq2_DT24_vs_DMSO_res1.csv")
res2.sort <- read.csv("DESeq2_DT48_vs_DMSO_res2.csv")


res2.out<-res2[complete.cases(res2),]
res2.out$minus_log10_adj.pval <- (-log10(res2.out$padj))
idx<-match(rownames(res2.out),counts.inf$Ensembl)
res2.out.inf<-cbind(counts.inf[idx,],data.frame(res2.out))
dim(res2.out.inf) #14762     11
res2.sort<-res2.out.inf[order(res2.out.inf$padj),]
rownames(res2.sort)<-NULL
write.csv(res2.sort,file="DESeq2_DT48_vs_DMSO_res2.csv")

########## MA plot 
pdf("plotMA_results.pdf", height=5, width=5)
plotMA(res1, ylim=c(-2,2),main="DESeq2_DT24_vs_DMSO_MAplot")
plotMA(res2, ylim=c(-2,2),main="DESeq2_DT48_vs_DMSO_MAplot")
dev.off()

#################


################Glimma
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Glimma")
library('Glimma')
Glimma::glimmaMA(res1.out,counts=counts.out)

glimmaMA(res2.out, dge=counts.out, status.colours=c("#3977db","#3d3f42","#db0d4e"))


#################
sig1<-res1.sort[!is.na(res1.sort$padj) & res1.sort$padj<0.05 & 
                  abs(res1.sort$log2FoldChange)>=1,]
head(sig1)
dim(sig1) #875 11

sig2<-res2.sort[!is.na(res2.sort$padj) & res2.sort$padj<0.05 & 
                  abs(res2.sort$log2FoldChange)>=1,]
head(sig2)
dim(sig2) #1205 11

selected1 <-rownames(sig1);selected1
selected2 <-rownames(sig2);selected2


library(gplots)
library(RColorBrewer)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDb")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")

library('SummarizedExperiment')

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) ## hmcol <- heat.colors

heatmap.2(log2(counts(res1,normalized=TRUE)[rownames(dds) %in% selected1,]),
          col=hmcol, scale="row",
          Rowv=TRUE, Colv=FALSE, 
          dendrogram="row",
          trace="none",
          margin=c(4,6), cexRow=0.5, cexCol=1, keysize=1)

counts(res1)
######### heatmap

library("pheatmap")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
select

sym<-counts.inf$HGNC[match(rownames(res2.sort)[1:20],counts.inf$Ensembl)]
ntd <- normTransform(dds)
rownames(ntd)<-sym

str(ntd)
vsd <- vst(dds, blind=FALSE)
rownames(vsd)<-sym

dfm <- as.data.frame(colData(dds)[,c("colnames.counts.","conditions")])
pheatmap(assay(vsd)[select,][complete.cases(rownames(assay(vsd)[select,])),], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=dfm)


pheatmap(assay(vsd)[res1.sort$HGNC[1:10],], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=dfm)

pheatmap(log1p(assay(ntd))[res1.sort$HGNC[1:10],], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=dfm)



head(counts)
dim(counts) #33158, 6
head(counts.out)
dim(counts.out) #33158     10
counts.out1 <- counts.out[!duplicated(counts.out$HGNC),] 
dim(counts.out1) #33100

counts<-counts.out1[,5:10]
head(counts)
dim(counts) #33100 6
head(counts.out1)
rownames(counts)<-counts.out1$HGNC
head(counts)


####################

### pre-processing of STN49 count matrix

counts_STN49 <-read.csv("gene_count.csv")
dim(counts_STN49) #58735 18

counts.out1 <- counts_STN49[!duplicated(counts_STN49$gene_name),] 
dim(counts.out1) #57169 18
counts<-counts.out1[,2:9]
dim(counts) #57169 8
head(counts)
rownames(counts)<-counts.out1$gene_name

class(counts)
counts[]<-lapply(counts,as.numeric)
counts[is.na(counts)]<-0

lcpm<-log2((counts+1)/colSums(counts+1)*1000000)
head(lcpm)

samples<-data.frame(colnames(counts))
samples<-data.frame(samples)
conditions<-colnames(counts)

samples<-data.frame(colnames(counts))
conditions<-c("CTRL24","CTRL24","BMi24","BMi24","CTRL48","CTRL48","BMi48","BMi48")
samples$conditions<-conditions
samples



#################

IFNg_genes <- c("GBP6","FGL2","EIF2AK2","IL2RB","STAT4","IRF1","NMI","UBE2L6","DHX58","USP18",
                "IFITM3","CASP1","PSMB9","NOD1","CIITA","SOD2","B2M","SRI","CXCL10","CXCL11",
                "SECTM1","IL6","SOCS3","CD74","DDX60","IFIT3","TRIM25","TRIM14","STAT1",
                "XAF1","MX1","OAS2","OAS3","SLC25A28","IFI35","IRF7","HLA-DQA1","RTP4","DDX58",
                "APOL6","PARP14","ZNFX1","SERPING1","RSAD2","HELZ2","CMPK2","STAT2","IFIH1",
                "TNFSF10","HERC6","C1R","HLA-DMA","C1S","HLA-A","HLA-B","IFI27","NFKBIA",
                "GBP4","HLA-DRB1","PARP12","ISG15","BST2","EPSTI1","PIM1","LGALS3BP","SAMHD1",
                "IFIT2")

temp <- lcpm[IFNg_genes,]
write.csv(temp,"human_IFNg_genes_list.csv")

pdf("heatmap_hs_IFN_genes.pdf", height=8, width=5)
heatmap.2(as.matrix(temp) ,main="IFNg_genes",
          scale="row", trace="none", key=TRUE,density.info=c("none"),
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          dendrogram="row", srtCol=45, cexCol =1, cexRow=0.5,
          margins = c(7, 8), Colv=FALSE,
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()



######## 


#################
Interest_genes1 <- c("NES","PROM1","GFAP","CLU","AQP4","CSPG4","PDGFRA","CNP","MBP",
                     "PLP1","PLLP","CLDN11","BTG2","ASPA","OLIG2","EGFR","HES6","ASCL1","DLL1","CD274","LGALS3")

temp <- lcpm[Interest_genes1,]

library(circlize)
library(ComplexHeatmap)

col_fun = colorRamp2(breaks=c(-2, 0, 2),colors=c("blue","white", "red"))
mat<-t(scale(t(as.matrix(temp))))

mat <-na.omit(mat)

# aGBM5
ave.mat<-cbind(( mat[,1]+mat[,2])/2,(mat[,3]+mat[,4])/2,(mat[,5]+mat[,6])/2)
# STN49 
ave.mat<-cbind(( mat[,1]+mat[,2])/2,(mat[,3]+mat[,4])/2,(mat[,5]+mat[,6])/2,(mat[,7]+mat[,8])/2)

pdf("heatmap_glial_genes_STN49_ave.pdf", height=10, width=20)
Heatmap(ave.mat, name = "Expression",  
        cluster_columns = TRUE,
        show_column_dend = TRUE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 1.5),
        column_gap = unit(0.5, "mm"),
        cluster_rows = TRUE,
        show_row_dend = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 8),
        column_title_rot = 45,
        show_column_names = TRUE,
        use_raster = FALSE,
        raster_quality = 20,
        width = unit(8, "cm"), height = unit(20, "cm"))
dev.off()



###################

IFNg_genes1 <-c("IFIT2","STAT2","SAMHD1","IFI27","IFITM2","IFITM3","MAPK3","FLNA",
               "IRF2","HLA-A","IFI6","IRF7","SOCS3","IRF5",
               "PSMB8","TRIM3","UBA7","IRF1","GBP3","OAS3","PLCG1","ADAR","XAF1","IFIT3",
               "ARIH1","IFNAR1","FLNB","MX2","NUP133","EIF4E2","MX1","RSAD2","USP18",
               "HERC5","IFIT1","CAMK2G","TRIM46","GBP1","BST2","IFITM1", "EIF4A1",
               "GBP2","HLA-F","OAS1", "STAT1","OAS2","JAK1","IP6K2","ISG15","EIF4A2","TPR",
               "EIF4E3","EIF4G3","UBE2L6","TRIM62","HLA-C","OASL","HLA-DQB1","IFNAR2","IRF4","TRIM2",
               "IFI35","IRF9","HLA-DRA","JAK2","HLA-DPA1","B2M","PIAS1","HLA-DRB5","IFNGR1","TRIM21",
               "TRIM22","TRIM14","CD44","SP100","SEC13","HLA-DRB1","NUP160","MID1","FCGR1A","PML",
               "NEDD4","TRIM68","TRIM17","NCAM1","CIITA","TRIM5","GBP5","GBP4","TRIM25","HLA-DQA1",
               "EIF2AK2","TRIM26","TRIM6","HLA-E","HLA-B")


Tcellactivation <- c("Cd4","Cd8a","Cd3","Foxp3","Tcf7","Cd69","Itgae","Itga1","Cd7","Il2ra",
                     "Tnf","Tbx21","Cd27","Cd28","Ifng","Cxcr6","Cd44","Gzmk","Gzmb","Prf1",
                     "Lck","Klrg1","Pdcd1","Tnfrsf4","Tnfrsf18","Tnfrsf9")
Tcellexhaustion <- c("Pdcd1","Lag3","Ctla4","Havcr2","Tigit","Btla","Eomes","Ccl3","Ccl4","Lgals3","Fas","Socs3")

Tcellactivation <-toupper(Tcellactivation)
Tcellexhaustion <-toupper(Tcellexhaustion)


MHCgenes <- c("HLA-DRB1","HLA-DRB3","HLA-DQA1","HLA-DQB1","HLA-DPB1","HLA-DMA","HLA-DMB","IFNGR1",
              "IFNGR2","IFNG","HLA-DOA","HLA-DOB","HLA-DQ","HLA-E","HLA-DQA1","HLA-DQB1","HLA-DQ",
              "HLA-DQA1","HLA-DQB1")


cytokines <- c("Ccrl2","Il6ra","Il6","Il1b","Ccl3","Ccl2","Ccl9","Cxcl2","Il1rn","Ccl6","Il17ra","Il22","Cxcl3","Il9", "Il17","Il23","Ifng","Tgfb",
               "Lif","Lifr","Mmp9")
cytokines <-toupper(cytokines)

write.csv(IFNg_genes1, "IFN_genes1_list.csv")


temp <- lcpm[Tcellexhaustion,]
temp <- na.omit(temp)

pdf("heatmap_hs_Tcell_exhaustion_genes.pdf", height=12, width=4)
heatmap.2(as.matrix(temp),main="Tcell_exhaustion_genes",
          scale="row", trace="none", key=TRUE,density.info=c("none"),
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          dendrogram="row", srtCol=45, cexCol =1, cexRow=0.6,
          margins = c(7, 8), Colv=FALSE,
          col=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

#########

pdf("Galectin3_expression_in_STN49.pdf", height=3, width=2)
stripchart(as.matrix(lcpm)[rownames(lcpm)%in%c("LGALS3"),]~samples$conditions,vertical=TRUE,las=2,cex.axis=1,pch=16,cex=1.3,
           col=mycol,method="jitter",ylab="expression",main=c("Galectin-3"))
dev.off()


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

####GO analysis 48hrs
human_genes_48hrs <- read.csv("DESeq2_DT48_vs_DMSO_res2.csv", header=TRUE)
head(human_genes_48hrs)
human_genes_48hrs_up <-subset(human_genes_48hrs,subset=human_genes_48hrs$log2FoldChange>0)
human_genes_48hrs_down <-subset(human_genes_48hrs,subset=human_genes_48hrs$log2FoldChange<0)

hs_genelist_up<-human_genes_48hrs_up$Entrez
hs_genelist_dn<-human_genes_48hrs_down$Entrez

hs_genelist_up<-list(X0=hs_genelist_up)
hs_genelist_dn<-list(X1=hs_genelist_dn)

lfc_up <-human_genes_48hrs_up$log2FoldChange
names(lfc_up) <-human_genes_48hrs_up$Entrez
lcf_up<-na.omit(lfc_up)
lcf_up<-sort(lcf_up,decreasing=TRUE)
head(lcf_up)

lfc_down <-human_genes_48hrs_down$log2FoldChange
names(lfc_down) <-human_genes_48hrs_down$Entrez
lcf_down<-na.omit(lfc_down)
lcf_down<-sort(lcf_down,decreasing=TRUE)

lcf <-c(lcf_up,lcf_down)
head(lcf)
tail(lcf)
lcf<-na.omit(lcf)

head(hs_genelist_up)
Combined_genelist<-cbind(hs_genelist_up,hs_genelist_dn)
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
ck_up<-compareCluster(geneClusters=Combined_genelist,fun="enrichGO", pvalueCutoff=0.05, OrgDb=org.Hs.eg.db, ont=ont[j])
head(as.data.frame(ck_up))
pdf("Dotplot_hs48_genes_GO_MF.pdf", width=5, height=10)
dotplot(ck_up, showCategory=20, font.size=8)
dev.off()



####
####

ont<-c("MF","CC","BP")
j=3
ego1 <- enrichGO(gene     = Combined_genelist$`1`,
                 OrgDb    = org.Hs.eg.db,
                 ont      = ont[j],
                 pAdjustMethod="BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = TRUE)


head(ego1)
#pdf(paste0("Heatplot_mmBP_up.pdf"),height=5,width=15)
#heatplot(ego1)
#dev.off()

pdf(paste0("Cnetplot_hsB48_P_up.pdf"),height=20,width=25)
#cnetplot(ego1,categorySize="p.adjust")
cnetplot(ego1,categorySize="p.adjust", colorEdge=TRUE, 
         foldChange=lcf_up,showCategory=5)
dev.off()
#######


d<-godata('org.Hs.eg.db', ont="BP")
ego2<-pairwise_termsim(ego1,method="Wang",semData=d)
ck1_up<-pairwise_termsim(ck_up,method="Wang",semData=d)

pdf("emapplot_hs48_BPa_up.pdf",height=10, width=10)
emapplot(ego2)
dev.off()

pdf("emapplot_hs48_BPb_up.pdf",height=10, width=10)
emapplot(ck1_up)
dev.off()

##############
ont<-c("MF","CC","BP")
j=3
ego3 <- enrichGO(gene     = Combined_genelist$`2`,
                 OrgDb    = org.Hs.eg.db,
                 ont      = ont[j],
                 pAdjustMethod="BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = TRUE)
#######
pdf(paste0("Cnetplot_hs48_BP_down1.pdf"),height=10,width=15)
#cnetplot(ego1,categorySize="p.adjust")
cnetplot(ego3,categorySize="p.adjust", colorEdge=TRUE, foldChange=lcf_down,showCategory=5)
dev.off()
#######


d<-godata('org.Hs.eg.db', ont="BP")
ego4<-pairwise_termsim(ego3,method="Wang",semData=d)
ck1_up<-pairwise_termsim(ck_up,method="Wang",semData=d)

pdf("emapplot_hs48_BPa_down.pdf",height=12, width=12)
emapplot(ego4)
dev.off()

pdf("emapplot_hs48_BPb_down.pdf",height=8, width=8)
emapplot(ck1_up)
dev.off()



#######
kk <- enrichKEGG(gene         = Combined_genelist$`1`,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

pdf("Dotplot_hs48_genes_KK_up.pdf", width=5, height=5)
dotplot(kk)
dev.off()

kkk<- setReadable(kk,OrgDb=org.Hs.eg.db,keyType="ENTREZID")

pdf("Cnetplot_hs48_genes_Endocytosis.pdf", width=10, height=10)
cnetplot(kkk, showCategory="Endocytosis", foldChange=lcf_up)
dev.off()

########
kk_down <- enrichKEGG(gene         = Combined_genelist$`2`,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)
head(kk)

pdf("Dotplot_hs48_genes_KK_down.pdf", width=5, height=5)
dotplot(kk_down)
dev.off()


####### REACTOME
#install.packages("reactome.db")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("reactome.db")

###### REACTOME up
library(ReactomePA)
x <- enrichPathway(gene=Combined_genelist$`1`,pvalueCutoff=0.05, readable=T, organism="human")
head(as.data.frame(x))

pdf("Barplot_hs48_genes_Reactome_up.pdf", width=4, height=4)
barplot(x, showCategory=10, font.size =6)
dev.off()

pdf("Dotplot_hs48_genes_Reactome_up.pdf", width=7, height=8)
dotplot(x, showCategory=10)
dev.off()


pdf("Cnetplot_hs48_genes_Reactome_up.pdf", width=25, height=25)
cnetplot(x, categorySize="pvalue", foldChange=lcf_up)
dev.off()

pdf("Treeplot_hs48_genes_Reactome_up.pdf", width=20, height=10)
treeplot(pairwise_termsim(x))
dev.off()

res <- compareCluster(Combined_genelist, fun="enrichPathway", organism="human")
pdf("Dotplot_hs48_genes_Reactome_Up_down.pdf", width=5, height=6)
dotplot(res)
dev.off()

########### REACTOME Down
x1 <- enrichPathway(gene=Combined_genelist$`2`,pvalueCutoff=0.05, readable=T, organism="human")
head(as.data.frame(x1))

pdf("Barplot_hs48_genes_Reactome_down.pdf", width=4, height=4)
barplot(x1, showCategory=10, font.size =6)
dev.off()

pdf("Dotplot_hs48_genes_Reactome_down.pdf", width=7, height=8)
dotplot(x1, showCategory=10)
dev.off()


pdf("Cnetplot_hs48_genes_Reactome_down.pdf", width=15, height=15)
cnetplot(x1, categorySize="pvalue", foldChange=lcf_down)
dev.off()

pdf("Treeplot_hs48_genes_Reactome_down.pdf", width=20, height=10)
treeplot(pairwise_termsim(x1))
dev.off()


############
y <- gsePathway(lcf, nPerm=10000,
                pvalueCutoff=0.2,
                pAdjustMethod="BH", verbose=FALSE, organism = "human")
yy<- setReadable(y,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
res <- as.data.frame(y)
head(y)

y2<-pairwise_termsim(y,method="JC")

pdf("emapplot_hs48_GSEA_REACTOME.pdf",height=10, width=20)
emapplot(y2, color="pvalue")
dev.off()

pdf("Dotplot_hs48_GSEA_REACTOME.pdf",height=5, width=6)
dotplot(y, showCategory=10)
dev.off()

write.csv(yy,file="hs48_GSEA_results_REACTOME.csv")



pdf("Ridgeplot_hs48_GSEA_reactome.pdf", width=6, height=18)
ridgeplot(yy)
dev.off()

pdf("Treeplot_hs48_GSEA_reactome.pdf", width=20, height=10)
treeplot(pairwise_termsim(y))
dev.off()


?ridgeplot
?cnetplot
?viewPathway



############

install.packages("msigdbr")
library(msigdbr)
msigdbr_show_species()
H <- msigdbr(species="Homo sapiens", category="H") %>%
  dplyr::select(gs_name,entrez_gene)

en1<-enricher(names(lcf_up), TERM2GENE = H)
en1<-setReadable(en1,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")

pdf("hs48_hallmark_results_a.pdf", height=2, width=4)
barplot(en1,showCategory = 20, font.size =5)
dotplot(en1,showCategory = 20, font.size=5,size=NULL)
dev.off()

?dotplot

pdf("hs48_hallmark_results_b.pdf", height=3, width=16)
treeplot(pairwise_termsim(en1), font.size=4)
dev.off()


################