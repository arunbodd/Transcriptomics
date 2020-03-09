library(scran)
library(scater)
library(DropletUtils)
library(SingleCellExperiment)
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)    
set.seed(1000)
library(Seurat)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
args <- commandArgs(trailingOnly = TRUE)

##Input arguments

DIR=args[1] ##Input directory
FILE1 = args[2] #Input File1, Sample1/subgroup
FILE2 = args[3] #Input FIle2, Sample2/subgroup
FILE3 = args[4] #Input File3, Sample3/subgroup
SN1 = args[5] #add cell ID (neg or pos)
SN2 = args[6] #add cell ID
SN3 = args[7] #add cell ID
subgroup = args[8] #give Sample1 or 2 or 3

setwd(DIR)
Sample_neg <- Read10X(data.dir = FILE1)
Sample_neg <- CreateSeuratObject(raw.data = Sample_neg, min.cells = 3, 
                                  min.genes = 200, project = paste0(SN1,"_",
                                                                    subgroup))

Sample_neg2 <- Read10X(data.dir = FILE2)
Sample_neg2 <- CreateSeuratObject(raw.data = Sample_neg2,min.cells = 3, 
                                  min.genes = 200, project = paste0(SN1,"_",
                                                                    subgroup))

Sample_neg3 <- Read10X(data.dir = FILE3)
Sample_neg3 <- CreateSeuratObject(raw.data = Sample_neg3,min.cells = 3, 
                                  min.genes = 200, project = paste0(SN1,"_",
                                                                    subgroup))

Sample.combined <- MergeSeurat(object1 = Sample_neg, object2 = Sample_neg2,
                                add.cell.id1 = SN1, 
                             add.cell.id2 = SN2,
                                project = paste0(subgroup))

Sample.combined <- AddSamples(Sample.combined, new.data = Sample_neg3, 
                              add.cell.id = SN3, project = paste0(subgroup))

table(Sample.combined@meta.data$orig.ident)
Sample_combined_sce = Convert(from = Sample.combined, to ="sce")

##Convert GeneSYMBOLS TO GENEIDs and remove NA
gene.ids <- mapIds(org.Hs.eg.db, keys=rownames(Sample_combined_sce),
                   keytype="SYMBOL", column="ENSEMBL")
keep <- !is.na(gene.ids) & !duplicated(gene.ids)
Sample_combined_sce <- Sample_combined_sce[keep,]
rownames(Sample_combined_sce) <- gene.ids[keep]
summary(keep)

##Get CHROMOSOME LOCATIONS
loc <- mapIds(EnsDb.Hsapiens.v86, keys=rownames(Sample_combined_sce),
              keytype="GENEID", 
              column="SEQNAME")
rowData(Sample_combined_sce)$Chr <- loc


##Pre-processing, removing lowlib, mitogenes
Sample_combined_sce2 <- calculateQCMetrics(Sample_combined_sce, compact=TRUE, 
                                      feature_controls=list(Mt=which(loc=="MT")))
QC <- Sample_combined_sce$scater_qc
lowlib <- isOutlier(Sample_combined_sce2$scater_qc$all$log10_total_counts, 
                    type="lower", nmads=3)
lowfeat <- isOutlier(Sample_combined_sce2$scater_qc$all$log10_total_features_by_counts, 
                     type="lower", nmads=3)
highmito <- isOutlier(Sample_combined_sce2$scater_qc$feature_control_Mt$pct_counts,
                      type="higher", nmads=3)
discard <- lowlib | highmito | lowfeat
Sample_combined_sce2 <- Sample_combined_sce2[,!discard]
summary(discard)

####We compute size factors for the endogenous genes and spike-in transcripts###
###and use them to compute log-normalized expression values####

set.seed(1000)
clusters <- quickCluster(Sample_combined_sce2, min.mean=0.1, method="igraph")
table(clusters)

Sample_combined_sce2 <- computeSumFactors(Sample_combined_sce2,
                                           min.mean=0.1, clusters=clusters)
saveRDS(file=paste0("FastMNNresults/",subgroup,".rds"), Sample_combined_sce2)
summary(sizeFactors(Sample_combined_sce2))
Sample_combined_sce2 <- normalize(Sample_combined_sce2)
fit <- trendVar(Sample_combined_sce2, parametric=TRUE,use.spikes=F) 
dec <- decomposeVar(Sample_combined_sce2, fit)
plot(dec$mean, dec$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)
dec.Sample_combined_sce2 <- dec
dec.Sample_combined_sce2$Symbol <- rowData(Sample_combined_sce2)$gene
dec.SSample_combined_sce2<- dec.Sample_combined_sce2[order
                                              (dec.Sample_combined_sce2$bio,
                                                         decreasing=TRUE),]

saveRDS(file=paste0("FastMNNresults/",subgroup,"_decomp",".rds"), dec.SSample_combined_sce2)
