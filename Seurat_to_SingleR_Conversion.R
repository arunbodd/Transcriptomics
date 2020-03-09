library(Seurat)
library(SingleR)
args <- commandArgs(trailingOnly = TRUE)

FILE1 = args[1]
subgroup = args[2]

load(file=("singleR_CSCE.rds"))
setwd("../")
singleR = readRDS(file = "singleR_CSCE.rds")
csce = readRDS(file = "FastMNNresults/CSCEobject.rds")
#Create a singleR object from Aligned data
singleR = CreateSinglerObject(logcounts(csce),project.name = "24h",
                              min.genes=0, technology = "10X",
                              species = "Human", normalize.gene.length = F, 
                              variable.genes = "de",
                              fine.tune = T, do.signatures = T, clusters=NULL, 
                              do.main.types = T,ref.list = HPCA,
                              reduce.file.size = T, numCores = 8)

save(singleR, file = paste0(subgroup,"_","SingleRobject.rds"))
#Add additional columns to metadata, new ident from Aligned sample, 
#tsne coordinates and resolution
singleR$meta.data$orig.ident = paste0(sce$Batch,"_",sce$Cluster)
singleR$meta.data$xy =  sce@reducedDims$MNN
singleR$meta.data$clusters = sce$Cluster


#Testing
singleR$meta.data$orig.ident = paste0(csce$Batch,"_",csce$Cluster)
singleR$meta.data$xy = csce@reducedDims$MNN
singleR$meta.data$clusters = csce$Cluster


SingleR.PlotTsne(singleR$singler[[1]]$SingleR.single.main,
                 singleR$meta.data$xy, do.label = T, do.letters = T, 
                 labels = singleR$meta.data$clusters, 
                 label.size = 2, dot.size = 3)

save(singleR, file = paste0("aggr24h","_","SingleRobject.rds"))

View(singleR$meta.data$xy)

#Draw TSNE plot for singleR object io
out = SingleR.PlotTsne(singleR$singler[[1]]$SingleR.single,
                       singleR$meta.data$xy, do.label = T, do.letters = T, 
                       labels = singleR$meta.data$clusters, 
                       label.size = 2, dot.size = 3)
pdf(paste0(subgroup,"_","SingleR_tSNEplot_clusters.pdf"), height=10, width=15)
out$p
dev.off()

pdf(paste0(subgroup,"_","SingleR_ClusterHeatmap.pdf"), height=10, width=15)
SingleR.DrawHeatmap(singleR$singler[[1]]$SingleR.single.main, top.n = 50, 
                    clusters = singleR$meta.data$clusters,
                    normalize = F,
                    order.by.clusters = T)
dev.off()

pdf(paste0(subgroup,"_","SingleR_groupHeatmap.pdf"), height = 10, width = 15)
SingleR.DrawHeatmap(singleR$singler[[1]]$SingleR.single.main, top.n = 10 , 
                    clusters = singleR$meta.data$orig.ident,
                    normalize = F,
                    order.by.clusters = T)
dev.off()
#Union_align = Seurat::StashIdent(object = Union_align,
#                                 save.name = "Cluster_ID")
#Color TSNE with cell subtypes
out = SingleR.PlotTsne(singleR$singler[[1]]$SingleR.single,
                       singleR$meta.data$xy, do.label = F, do.letters = F, 
                       labels = singleR$singler[[1]]$SingleR.single$labels, 
                       label.size = 5, dot.size = 2)
pdf(paste0(subgroup,"_","tSNE-cellsubtypes.pdf"), height = 10, width = 15)
out$p
dev.off()

#Output the COUNTS per cluster for each group
html(paste0(subgroup,"_","Clustercounts.html"), height=10, width =15)
x = knitr::kable(table(singleR$singler[[1]]$SingleR.single$labels,
                       singleR$meta.data$orig.ident),"latex")
x1=kable_styling(x, font_size = 6, bootstrap_options = "striped", full_width = F)
as_image(x1,file = paste0(subgroup,"_","Clustercounts.jpeg"))
