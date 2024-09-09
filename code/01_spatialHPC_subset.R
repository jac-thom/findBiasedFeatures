library(SpatialExperiment)
library(dplyr)

load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/06_clustering/PRECAST/spe_precast_HE_domain.rda')

fix_order = distinct(as.data.frame(colData(spe)), slide, array, brnum, sample_id, position, sex) %>% 
  arrange(slide, array)
sub4 = fix_order$sample_id[c(14,16,
                             20,21)]
spe_sub4 = spe[,spe$sample_id %in% sub4]

#load and reformat Erik's nnSVG results
load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/nnSVG/nnSVG_outs_HE_only.rda")
res_df = tidyr::pivot_longer(tibble::rownames_to_column(as.data.frame(res_ranks), var="gene_id"), 
                             colnames(res_ranks), names_to="sample_id", values_to="rank", values_drop_na=T)
#filter to only the top 2k sig features in the 4 samples we're using
res_df2 = filter(res_df, sample_id %in% c("V11L05-333_B1","V11L05-333_D1","V11L05-335_D1","V11L05-336_A1"),
       rank<=2000)
#further filter to only features that are in the top 2k of >1 sample
svgs = group_by(res_df2, gene_id) %>% tally() %>% filter(n>1) 
nrow(svgs)

spe_sub4 = spe_sub4[svgs$gene_id,]
dim(spe_sub4)

write.csv(counts(spe_sub4)[1:1000,], "/users/jthompso/BiasedFeatures_hpc/final/spe-hpc_sub4_svgs-only_counts-1.csv", row.names=T)
write.csv(counts(spe_sub4)[1001:2082,], "/users/jthompso/BiasedFeatures_hpc/final/spe-hpc_sub4_svgs-only_counts-2.csv", row.names=T)
write.csv(colData(spe_sub4)[,c("sample_id","slide","brnum","array","sex","in_tissue","array_row","array_col","sum_umi","sum_gene","expr_chrM","expr_chrM_ratio")],
          "/users/jthompso/BiasedFeatures_hpc/final/spe-hpc_sub4_svgs-only_colData.csv", row.names=T)
write.csv(rowData(spe_sub4), "/users/jthompso/BiasedFeatures_hpc/final/spe-hpc_sub4_svgs-only_rowData.csv", row.names=T)
