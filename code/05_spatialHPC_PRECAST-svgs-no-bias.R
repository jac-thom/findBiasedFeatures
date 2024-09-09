suppressPackageStartupMessages({
  library(Seurat)
  library(PRECAST)
  library(dplyr)
  library(here)
})
set.seed(123)

#build spe
counts = rbind(read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_counts-1.csv", row.names=1),
               read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_counts-2.csv", row.names=1))
colnames(counts) = gsub("\\.","-",colnames(counts))
cdata = read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_colData.csv", row.names=1)
rdata = read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_rowData.csv", row.names=1)

spe = SpatialExperiment(
  assay = list("counts"=counts), 
  colData = cdata, rowData = rdata,
  spatialCoordsNames = c("array_row", "array_col"))

#reformat spe to seurat list
l2 = unique(spe$sample_id)
names(l2) = l2
l2 = lapply(l2, function(x) spe[,colData(spe)$sample_id==x])

srt.sets = lapply(l2, function(x) {
  colnames(counts(x)) <- rownames(colData(x))
  colData(x)$col <- x$array_col
  colData(x)$row <- x$array_row
  count <- counts(x)
  a1 <- CreateAssayObject(count, assay = "RNA", min.features = 0, min.cells = 0)
  CreateSeuratObject(a1, meta.data = as.data.frame(colData(x)))
})


#remove biased genes
bias = readRDS("findBiasedFeatures/processed-data/hpc_biased-features.rda")
svgs_filt = setdiff(rownames(spe), bias)
length(svgs_filt) #should be 2067

#run precast
preobj <- CreatePRECASTObject(seuList = srt.sets,
                              customGenelist=svgs_filt,
                              premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj <- AddAdjList(preobj, platform = "Visium")
PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj <- PRECAST(PRECASTObj, K=7)

#consolidate/ reformat results
PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

write.csv(seuInt@meta.data, "findBiasedFeatures/processed-data/seuInt-hpc_k-7_svgs-no-bias_metadata.csv", row.names=T)