library(SpatialExperiment)
library(scater)
library(ggspavis)

#build spe
counts = rbind(read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_counts-1.csv", row.names=1),
               read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_counts-2.csv", row.names=1))
colnames(counts) = gsub("\\.","-",colnames(counts))
cdata = read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_colData.csv", row.names=1)
rdata = read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_rowData.csv", row.names=1)
spdata = read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_spatialCoords.csv", row.names=1)

spe = SpatialExperiment(
  assay = list("counts"=counts), 
  colData = cdata, rowData = rdata,
  spatialCoords = as.matrix(spdata))

#add PRECAST results to colData
clusters1 = read.csv("findBiasedFeatures/processed-data/seuInt-hpc_k-7_svgs_metadata.csv", row.names=1)
identical(rownames(clusters1), rownames(colData(spe)))
spe$precast_k7 = clusters1$cluster
spe$precast_k7_ordered = factor(spe$precast_k7, levels=c(7,2,3,1,6,5,4), 
                                labels=c("WM","WM (2)","SR/SL","CA1","CA3","DG GCL","DG ML"))

clusters2 = read.csv("findBiasedFeatures/processed-data/seuInt-hpc_k-7_svgs-no-bias_metadata.csv", row.names=1)
identical(rownames(clusters2), rownames(colData(spe)))
spe$precast_k7_nobias = clusters2$cluster
spe$precast_k7_nobias_ordered= factor(spe$precast_k7_nobias, levels=c(1,2,7,5,6,4,3),
                                      labels=c("WM","SR/SL","CA1","CA1 (2)","CA3","DG GCL","DG ML"))

#heatmap to justify cluster annotations
spe <- scuttle::logNormCounts(spe)
markers = c("MBP","GFAP","SPARCL1","FIBCD1","COL5A2","KCNQ5","CARTPT","PCDH8","CALB1")

plotGroupedHeatmap(spe, features = markers, swap_rownames="gene_name", 
                   group="precast_k7_ordered",
                   scale=TRUE, center=TRUE, 
                   cluster_rows=FALSE, cluster_cols=FALSE)
plotGroupedHeatmap(spe, features = markers, swap_rownames="gene_name",
                   group="precast_k7_nobias_ordered",
                   scale=TRUE, center=TRUE, 
                   cluster_rows=FALSE, cluster_cols=FALSE)

#plot cluster results
l2 = unique(spe$sample_id)
names(l2) = l2
l2 = lapply(l2, function(x) spe[,colData(spe)$sample_id==x])

col.pal1 = c("#1f77b4FF","#aec7e8FF","#ffbb78FF","#2ca02cFF","#ff7f0eFF","#d62728FF","#ff9896FF")
col.pal2 = c("#1f77b4FF","#ffbb78FF","#2ca02cFF","#98df8aFF","#ff7f0eFF","#d62728FF","#ff9896FF")

c1 <- lapply(seq_along(l2), function(x) {
  plotSpots(l2[[x]], annotate="precast_k7_ordered", point_size=.3)+
    labs(color="clus")+
    scale_color_manual(values=col.pal1)+
    theme(plot.title=element_text(size=8))
})
PRECAST::drawFigs(c1, layout.dim = c(1, 4), common.legend = TRUE, legend.position = "right", align = "h")

c2 <- lapply(seq_along(l2), function(x) {
  plotSpots(l2[[x]], annotate="precast_k7_nobias_ordered", point_size=.3)+
    labs(color="clus")+
    scale_color_manual(values=col.pal2)+
    theme(plot.title=element_text(size=8))
})
PRECAST::drawFigs(c2, layout.dim = c(1, 4), common.legend = TRUE, legend.position = "right", align = "h")
