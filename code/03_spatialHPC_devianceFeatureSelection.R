library(SpatialExperiment)
library(dplyr)
library(scry)

counts = rbind(read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_counts-1.csv", row.names=1),
               read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_counts-2.csv", row.names=1))
colnames(counts) = gsub("\\.","-",colnames(counts))
cdata = read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_colData.csv", row.names=1)
rdata = read.csv("findBiasedFeatures/processed-data/spe-hpc_sub4_svgs-only_rowData.csv", row.names=1)

spe = SpatialExperiment(
  assay = list("counts"=counts), 
  colData = cdata, rowData = rdata,
  spatialCoordsNames = c("array_row", "array_col"))

#if need to visualize gene expression 
#spe <- scuttle::logNormCounts(spe)

default <- devianceFeatureSelection(spe, fam="binomial", batch=NULL)
df = cbind.data.frame("gene"=rownames(default),"gene_name"=rowData(default)$gene_name,
                      "dev"= rowData(default)$binomial_deviance,
                      "rank"=(nrow(default)+1)-rank(rowData(default)$binomial_deviance))

batch <- devianceFeatureSelection(spe, fam="binomial", batch=as.factor(spe$sample_id))

df = left_join(df, cbind.data.frame("gene"=rownames(rowData(batch)), 
                                    "dev"=rowData(batch)$binomial_deviance,
                                    "rank"=(nrow(batch)+1)-rank(rowData(batch)$binomial_deviance)),
               by="gene", suffix=c("_default","_sample"))

write.csv(df, "findBiasedFeatures/processed-data/hpc_bindev_default-sample_svgs-only.csv", row.names=FALSE)
